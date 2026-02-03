import re
import argparse
import sys
import os
import shutil
from datetime import datetime


def normalize_name(name):
    return name.replace("-", "__")


def create_backup(header_file):
    """Create a backup of the header file with a timestamp."""
    base_name, ext = os.path.splitext(header_file)  # Split base name and extension
    timestamp = int(datetime.now().timestamp())  # Current time in seconds
    backup_file = f"{base_name}_{timestamp}{ext}"  # Correctly format the backup file name
    shutil.copy2(header_file, backup_file)
    print(f"Backup created: {backup_file}")
    return backup_file


def categorize_reactions(species_file, reactions_list):
    """Parse the species file to extract reactions and their formulas."""
    reactions = {}
    with open(species_file, 'r') as sf:
        lines = sf.readlines()

    for line in lines:
        # Parse reaction formulas
        match = re.match(r'chemistry\.([\w\-]+)\.reaction\s*=\s*(.+)', line)
        if match:
            reaction_name = match.group(1)
            formula = match.group(2)
            reactions[reaction_name] = formula.strip()

    """Categorize reactions into lagrangian and eulerian based on their formulas."""
    lagrangian = []
    eulerian = []

    for reaction in reactions_list:
        formula = reactions.get(reaction, "")
        if "(s)" in formula:
            lagrangian.append((reaction, formula))
        else:
            eulerian.append((reaction, formula))

    return lagrangian, eulerian


def insert_optimized_usr_rates(header_file, state):
    """
    Insert '#define OPTIMIZED_USR_RATES' if not already present,
    below '#define MFIX_USR_REACTIONS_RATES_K_H'.
    
    The 'state' dictionary will store whether the file was already optimized.
    """
    with open(header_file, 'r') as file:
        lines = file.readlines()

    # Check if the line is already present
    if any(line.strip() == "#define OPTIMIZED_USR_RATES" for line in lines):
        state['already_optimized'] = 1
        return  # Exit if the line is already present

    # Locate the target line and insert below it
    for i, line in enumerate(lines):
        if line.strip() == "#define MFIX_USR_REACTIONS_RATES_K_H":
            lines.insert(i + 1, "#define OPTIMIZED_USR_RATES\n")
            state['already_optimized'] = 0
            break
    else:
        print("Error: #define MFIX_USR_REACTIONS_RATES_K_H not found. No changes made.")
        sys.exit(1)

    # Write the updated content back to the file
    with open(header_file, 'w') as file:
        file.writelines(lines)


def replace_species_idx(header_file, valid_species, species_type):
    """Replace fluid.species_idx("species_name") or solids.species_idx("species_name") 
    with fluid_idx::species_name or solids_idx::species_name in the header file."""
    if species_type not in ["fluid", "solids"]:
        raise ValueError(f"Invalid species_type: {species_type}. Must be 'fluid' or 'solids'.")

    with open(header_file, 'r') as file:
        content = file.read()

    # Construct patterns and replacements based on species_type
    pattern = rf'{species_type}\.species_idx\("([\w\-]+)"\)'
    replacement_prefix = f"{species_type}_idx::"

    # Function to perform conditional replacement
    def replacement_idx_function(match):
        species_name = match.group(1)
        if species_name in valid_species:
            return f"{replacement_prefix}{normalize_name(species_name)}"
        else:
            print(f"Error: '{species_name}' is not a valid {species_type} species name.")
            sys.exit(1)

    # Perform the replacement
    updated_content = re.sub(pattern, replacement_idx_function, content)

    # Construct patterns and replacements based on species_type
    pattern = rf'{species_type}\.(\w+)\("([\w\-]+)"\)'  # Match species_type.FUNCTION("species_name")
    replacement_prefix = f"{species_type}_idx::"

    # Function to perform conditional replacement
    def replacement_name_function(match):
        function_name = match.group(1)  # Extract FUNCTION
        species_name = match.group(2)  # Extract species_name
        if species_name in valid_species:
            return f"{species_type}.{function_name}({replacement_prefix}{normalize_name(species_name)})"
        else:
            print(f"Error: '{species_name}' is not a valid {species_type} species name.")
            sys.exit(1)

    # Perform the replacement
    final_content = re.sub(pattern, replacement_name_function, updated_content)

    # Write the updated content back to the header file
    with open(header_file, 'w') as file:
        file.write(final_content)


def replace_reaction_idx(header_file, valid_reactions):
    """Replace reactions.reaction_idx("reaction_name")
    with reaction_idx::reaction_name"""

    with open(header_file, 'r') as file:
        content = file.read()

    # Construct patterns and replacements based on species_type
    pattern = rf'reactions.reaction_idx\("([\w\-]+)"\)'
    replacement_prefix = f"reactions_idx::"

    # Function to perform conditional replacement
    def replacement_idx_function(match):
        reaction_name = match.group(1)
        if reaction_name in valid_reactions:
            return f"{replacement_prefix}{normalize_name(reaction_name)}"
        else:
            print(f"Error: '{reaction_name}' is not a valid reaction.")
            sys.exit(1)

    # Perform the replacement
    updated_content = re.sub(pattern, replacement_idx_function, content)

    # Construct patterns and replacements based on species_type
    pattern = rf'reactions\.(\w+)\("([\w\-]+)"\)'  # Match reactions.FUNCTION("species_name")
    replacement_prefix = f"reactions_idx::"

    # Function to perform conditional replacement
    def replacement_name_function(match):
        function_name = match.group(1)  # Extract FUNCTION
        reaction_name = match.group(2)  # Extract reaction_name
        if reaction_name in valid_reactions:
            return f"reactions.{function_name}({replacement_prefix}{normalize_name(reaction_name)})"
        else:
            print(f"Error: '{reaction_name}' is not a valid reactions name.")
            sys.exit(1)

    # Perform the replacement
    final_content = re.sub(pattern, replacement_name_function, updated_content)

    # Write the updated content back to the header file
    with open(header_file, 'w') as file:
        file.write(final_content)



def add_struct_enums(header_file, fluid_species_list, solids_species_list, lagrangian, eulerian):
    # Read the header file content
    with open(header_file, 'r') as hf:
        header_lines = hf.readlines()

    # Process the content to insert the structs at the beginning of the
    # LagrangianReactionRates and EulerianReactionRates classes
    processed_lines = []
    in_lag_class = False
    in_eul_class = False
    for line in header_lines:
        if re.match(r'\bclass\b\s+LagrangianReactionRates', line):
            in_lag_class = True
            processed_lines.append(line)
        elif in_lag_class and "{" in line:  # Opening brace of the class
            processed_lines.append(line)
            processed_lines.append("  public:\n")

            # Generate and insert the struct content
            if fluid_species_list:
                fluid_struct_content = "    struct fluid_idx { enum {\n"
                for idx, species in enumerate(fluid_species_list):
                    if idx != len(fluid_species_list) - 1:  # Execute only if idx is not the last index
                        fluid_struct_content += f"      {normalize_name(species)} = {idx},\n"
                    else:
                        fluid_struct_content += f"      {normalize_name(species)} = {idx}\n"
                fluid_struct_content += "    };};\n\n"
                processed_lines.append(fluid_struct_content)

            # Generate and insert the struct content
            if solids_species_list:
                solids_struct_content = "    struct solids_idx { enum {\n"
                for idx, species in enumerate(solids_species_list):
                    if idx != len(solids_species_list) - 1:  # Execute only if idx is not the last index
                        solids_struct_content += f"      {normalize_name(species)} = {idx},\n"
                    else:
                        solids_struct_content += f"      {normalize_name(species)} = {idx}\n"
                solids_struct_content += "    };};\n\n"
                processed_lines.append(solids_struct_content)

            # Create struct for lagrangian reactions
            if lagrangian:
                reactions_struct_content = "    struct reactions_idx { enum {\n"
                for idx, (name, _) in enumerate(lagrangian):
                    if idx != len(lagrangian) - 1:  # Execute only if idx is not the last index
                        reactions_struct_content += f"      {normalize_name(name)} = {idx},\n"
                    else:
                        reactions_struct_content += f"      {normalize_name(name)} = {idx}\n"
                reactions_struct_content += "    };};\n\n"
                processed_lines.append(reactions_struct_content)

            in_lag_class = False  # No need to track class body further for this operation
        elif re.match(r'\bclass\b\s+EulerianReactionRates', line):
            in_eul_class = True
            processed_lines.append(line)
        elif in_eul_class and "{" in line:  # Opening brace of the class
            processed_lines.append(line)
            processed_lines.append("  public:\n")

            # Generate and insert the struct content
            if fluid_species_list:
                fluid_struct_content = "    struct fluid_idx { enum {\n"
                for idx, species in enumerate(fluid_species_list):
                    if idx != len(fluid_species_list) - 1:  # Execute only if idx is not the last index
                        fluid_struct_content += f"      {normalize_name(species)} = {idx},\n"
                    else:
                        fluid_struct_content += f"      {normalize_name(species)} = {idx}\n"
                fluid_struct_content += "    };};\n\n"
                processed_lines.append(fluid_struct_content)

            # Create struct for eulerian reactions
            if eulerian:
                reactions_struct_content = "    struct reactions_idx { enum {\n"
                for idx, (name, _) in enumerate(eulerian):
                    if idx != len(eulerian) - 1:  # Execute only if idx is not the last index
                        reactions_struct_content += f"      {normalize_name(name)} = {idx},\n"
                    else:
                        reactions_struct_content += f"      {normalize_name(name)} = {idx}\n"
                reactions_struct_content += "    };};\n\n"
                processed_lines.append(reactions_struct_content)

            in_eul_class = False  # No need to track class body further for this operation
        else:
            processed_lines.append(line)

    # Write the updated content back to the header file
    with open(header_file, 'w') as file:
        file.write("".join(processed_lines))



def add_constexpr_arrays(header_file, fluid_species_list, solids_species_list, lagrangian_list, eulerian_list):
    # Read the header file content
    with open(header_file, 'r') as hf:
        header_lines = hf.readlines()

    # Prepare the species array declaration
    fluid_species_count = len(fluid_species_list)
    solids_species_count = len(solids_species_list)
    lagrangian_count = len(lagrangian_list)
    eulerian_count = len(eulerian_list)

    fluid_array_declaration = ""
    if fluid_species_list:
        fluid_array_declaration = (
            "\nnamespace rates_aux { namespace fluid_aux {\n"
            + f"static amrex::Array<std::string, {fluid_species_count}> species_names = {{\n"
            + ",\n".join(f'    "{normalize_name(species)}"' for species in fluid_species_list)
            + "\n};\n}} // end namespace rates_aux::fluid_aux\n\n"
        )

    solids_array_declaration = ""
    if solids_species_list:
        solids_array_declaration = (
            "\nnamespace rates_aux { namespace solids_aux {\n"
            + f"static amrex::Array<std::string, {solids_species_count}> species_names = {{\n"
            + ",\n".join(f'    "{normalize_name(species)}"' for species in solids_species_list)
            + "\n};\n}} // end namespace rates_aux::solids_aux\n\n"
        )

    lagrangian_array_declaration = ""
    if lagrangian_list:
        lagrangian_array_declaration = (
            "\nnamespace rates_aux { namespace lagrangian_aux {\n"
            + f"static amrex::Array<std::string, {lagrangian_count}> reactions_names = {{\n"
            + ",\n".join(f'    "{normalize_name(reaction)}"' for reaction in lagrangian_list)
            + "\n};\n}} // end namespace rates_aux::lagrangian_aux\n\n"
        )

    eulerian_array_declaration = ""
    if eulerian_list:
        eulerian_array_declaration = (
            "\nnamespace rates_aux { namespace eulerian_aux {\n"
            + f"static amrex::Array<std::string, {eulerian_count}> reactions_names = {{\n"
            + ",\n".join(f'    "{normalize_name(reaction)}"' for reaction in eulerian_list)
            + "\n};\n}} // end namespace rates_aux::eulerian_aux\n\n"
        )

    # Insert the species array before the **last** #endif directive
    final_content = []
    last_endif_index = None
    for idx, line in enumerate(header_lines):
        if line.strip() == "#endif":
            last_endif_index = idx  # Track the index of the last #endif

    if last_endif_index is not None:
        final_content = (
            header_lines[:last_endif_index]  # All lines before the last #endif
            + [fluid_array_declaration]         # Insert the array declaration
            + [solids_array_declaration]        # Insert the array declaration
            + [lagrangian_array_declaration]    # Insert the array declaration
            + [eulerian_array_declaration]      # Insert the array declaration
            + [header_lines[last_endif_index]]  # Add the last #endif back
        )
    else:
        final_content = header_lines  # No #endif directive found, just use the lines as is

    with open(header_file, 'w') as of:
        of.write("".join(final_content))


def write_all_test_indexes(header_file, fluid_species_list, solids_species_list, lagrangian_list, eulerian_list):
    with open(header_file, 'r') as hf:
        header_lines = hf.readlines()

    fluid_count = len(fluid_species_list)
    solids_count = len(solids_species_list)
    lagrangian_count = len(lagrangian_list)
    eulerian_count = len(eulerian_list)

    first_part = (
        "\nnamespace rates_aux {\n"
        + "\ninline std::string normalize_name(const std::string& name) {\n"
        + "    std::string result = name;\n"
        + "    size_t pos = 0;\n"
        + "    while ((pos = result.find('-', pos)) != std::string::npos) {\n"
        + "        result.replace(pos, 1, \"__\");\n"
        + "        pos += 2;\n"
        + "    }\n"
        + "    return result;\n"
        + "}\n"
    )

    middle_parts = []

    # Fluid test function
    if fluid_species_list:
        middle_parts.append("\ntemplate <class Fluid>\nstatic int test_fluid_indexes (const Fluid& fluid)\n{\n")
        middle_parts.append("  auto const& fluid_parms = fluid.template parameters<amrex::RunOn::Host>();\n")
        middle_parts.append("  auto const& species_names = fluid.species_names();\n\n")
        middle_parts.append("  if (species_names.size() != fluid_aux::species_names.size()) return 1;\n\n")
        for i in range(fluid_count):
            s = fluid_species_list[i]
            sn = f"species_names[{i}]"
            middle_parts.append("  {\n")
            middle_parts.append(f"    const std::string& species_name = {sn};\n")
            middle_parts.append(f"    if (normalize_name(species_name) != fluid_aux::species_names[{i}]) return 1;\n")
            middle_parts.append(f"    if (fluid_parms.species_idx(species_name.c_str()) != LagrangianReactionRates::fluid_idx::{s.replace('-', '__')}) return 1;\n")
            middle_parts.append(f"    if (fluid_parms.species_idx(species_name.c_str()) != EulerianReactionRates::fluid_idx::{s.replace('-', '__')}) return 1;\n")
            middle_parts.append("  }\n")
        middle_parts.append("  return 0;\n}\n")
    else:
        middle_parts.append("\ntemplate <class Fluid>\nstatic int test_fluid_indexes (const Fluid&)\n{\n")
        middle_parts.append("  return 0;\n}\n")

    # Solids test function
    if solids_species_list:
        middle_parts.append("\ntemplate <class Solids>\nstatic int test_solids_indexes (const Solids& solids)\n{\n")
        middle_parts.append("  auto& solids_parms = solids.template parameters<amrex::RunOn::Host>();\n")
        middle_parts.append("  auto const& species_names = solids.species_names();\n\n")
        middle_parts.append("  if (species_names.size() != solids_aux::species_names.size()) return 1;\n\n")
        for i in range(solids_count):
            s = solids_species_list[i]
            sn = f"species_names[{i}]"
            middle_parts.append("  {\n")
            middle_parts.append(f"    const std::string& species_name = {sn};\n")
            middle_parts.append(f"    if (normalize_name(species_name) != solids_aux::species_names[{i}]) return 1;\n")
            middle_parts.append(f"    if (solids_parms.species_idx(species_name.c_str()) != LagrangianReactionRates::solids_idx::{s.replace('-', '__')}) return 1;\n")
            middle_parts.append("  }\n")
        middle_parts.append("  return 0;\n}\n")
    else:
        middle_parts.append("\ntemplate <class Solids>\nstatic int test_solids_indexes (const Solids&)\n{\n")
        middle_parts.append("  return 0;\n}\n")

    # Lagrangian reactions test function
    if lagrangian_list:
        middle_parts.append("\ntemplate <class Reactions>\nstatic int test_lagrangian_indexes (const Reactions& reactions)\n{\n")
        middle_parts.append("  auto& lagrangian_parms = reactions.template lagrangian_parameters<amrex::RunOn::Host>();\n\n")
        middle_parts.append("  const int nb_lagrangian = reactions.lagrangian_nreactions();\n")
        middle_parts.append("  if (nb_lagrangian != lagrangian_aux::reactions_names.size()) return 1;\n\n")
        for i in range(lagrangian_count):
            r = lagrangian_list[i]
            rn = f"reactions.get_lagrangian({i})->get_name()"
            middle_parts.append("  {\n")
            middle_parts.append(f"    const std::string& reaction_name = {rn};\n")
            middle_parts.append(f"    if (normalize_name(reaction_name) != lagrangian_aux::reactions_names[{i}]) return 1;\n")
            middle_parts.append(f"    if (lagrangian_parms.reaction_idx(reaction_name.c_str()) != LagrangianReactionRates::reactions_idx::{r.replace('-', '__')}) return 1;\n")
            middle_parts.append("  }\n")
        middle_parts.append("  return 0;\n}\n")
    else:
        middle_parts.append("\ntemplate <class Reactions>\nstatic int test_lagrangian_indexes (const Reactions&)\n{\n")
        middle_parts.append("  return 0;\n}\n")


    # Eulerian reactions test function
    if eulerian_list:
        middle_parts.append("\ntemplate <class Reactions>\nstatic int test_eulerian_indexes (const Reactions& reactions)\n{\n")
        middle_parts.append("  auto& eulerian_parms = reactions.template eulerian_parameters<amrex::RunOn::Host>();\n\n")
        middle_parts.append("  const int nb_eulerian = reactions.eulerian_nreactions();\n")
        middle_parts.append("  if (nb_eulerian != eulerian_aux::reactions_names.size()) return 1;\n\n")
        for i in range(eulerian_count):
            r = eulerian_list[i]
            rn = f"reactions.get_eulerian({i})->get_name()"
            middle_parts.append("  {\n")
            middle_parts.append(f"    std::string reaction_name = {rn};\n")
            middle_parts.append(f"    if (normalize_name(reaction_name) != eulerian_aux::reactions_names[{i}]) return 1;\n")
            middle_parts.append(f"    if (eulerian_parms.reaction_idx(reaction_name.c_str()) != EulerianReactionRates::reactions_idx::{r.replace('-', '__')}) return 1;\n")
            middle_parts.append("  }\n")
        middle_parts.append("  return 0;\n}\n")
    else:
        middle_parts.append("\ntemplate <class Reactions>\nstatic int test_eulerian_indexes (const Reactions&)\n{\n")
        middle_parts.append("  return 0;\n}\n")


    last_part = "\n} // end namespace rates_aux\n\n"

    final_content = []
    last_endif_index = None
    for idx, line in enumerate(header_lines):
        if line.strip() == "#endif":
            last_endif_index = idx

    if last_endif_index is not None:
        final_content = (
            header_lines[:last_endif_index]
            + [first_part]
            + middle_parts
            + [last_part]
            + [header_lines[last_endif_index]]
        )
    else:
        final_content = header_lines

    with open(header_file, 'w') as of:
        of.write("".join(final_content))


def update_header_file(species_file, header_file):
    """Update the header file based on the species definitions."""
    temp_file = f"{header_file}.temp"

    try:
        # Create a temporary copy of the header file
        shutil.copy2(header_file, temp_file)
        
        # State to track optimization status
        state = {'already_optimized': 0}

        # Insert the optimized rates definition and update state
        insert_optimized_usr_rates(temp_file, state)    

        # Read species from the input file
        with open(species_file, 'r') as sf:
            content = sf.read()

        # Remove line breaks after the continuation character '\'
        content = re.sub(r'\\\s*\n', ' ', content)

        # Parse the species names
        fluid_match = re.search(r'fluid\.species\s*=\s*(.+)', content)
        solids_match = re.search(r'solids\.species\s*=\s*(.+)', content)
        if not (fluid_match or solids_match):
            raise ValueError("Species information not found in the file.")

        # Parse the species names
        reactions_match = re.search(r'chemistry\.solve\s*=\s*(.+)', content)
        if not (reactions_match):
            raise ValueError("Reactions information not found in the file.")

        fluid_species_list = fluid_match.group(1).split() if fluid_match else []
        solids_species_list = solids_match.group(1).split() if solids_match else []
        reactions_list = reactions_match.group(1).split()

        # Categorize reactions
        lagrangian, eulerian = categorize_reactions(species_file, reactions_list)

        lagrangian_list = [inner_list[0] for inner_list in lagrangian]
        eulerian_list = [inner_list[0] for inner_list in eulerian]

        """Update the header file based on the species definitions."""
        replace_species_idx(temp_file, fluid_species_list, "fluid")
        replace_species_idx(temp_file, solids_species_list, "solids")
        replace_reaction_idx(temp_file, reactions_list)

        if not (state['already_optimized']):
            add_struct_enums(temp_file, fluid_species_list, solids_species_list, lagrangian, eulerian)
            add_constexpr_arrays(temp_file, fluid_species_list, solids_species_list, lagrangian_list, eulerian_list)
            write_all_test_indexes(temp_file, fluid_species_list, solids_species_list, lagrangian_list, eulerian_list)

        # Create a backup of the original file before finalizing changes
        create_backup(header_file)

        # Replace the original header file with the optimized file
        shutil.move(temp_file, header_file)
        print(f"Header file updated and saved to {header_file}")

    except SystemExit as e:
        print(f"SystemExit occurred with code {e.code}. Cleaning up...")
        if os.path.exists(temp_file):
            os.remove(temp_file)  # Clean up the temporary file
        sys.exit(e.code)  # Re-raise SystemExit to propagate the exit code

    except Exception as e:
        print(f"Error occurred: {e}. Original file remains unchanged.")
        if os.path.exists(temp_file):
            os.remove(temp_file)  # Clean up the temporary file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update a C++ header file based on species definitions.")
    parser.add_argument("species_file", help="Path to the file containing species definitions.")
    parser.add_argument("header_file", help="Path to the header file to be updated.")

    args = parser.parse_args()

    species_file = args.species_file
    header_file = args.header_file

    # Ensure the species and header files exist
    if not os.path.exists(species_file):
        print(f"Error: The file '{species_file}' does not exist.")
        sys.exit(1)

    if not os.path.exists(header_file):
        print(f"Error: The file '{header_file}' does not exist.")
        sys.exit(1)

    # Update the header file in place
    print(f"Updating header file: {header_file}")
    update_header_file(species_file, header_file)
