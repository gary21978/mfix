#ifndef _MFIX_EB_NODE_CONNECTIVITY_K_H_
#define _MFIX_EB_NODE_CONNECTIVITY_K_H_

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 0                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node0_EB(const int covered,
                        const amrex::IntVect Node0, const int /*Node0_isNotConnected*/,
                        const amrex::IntVect Node1, const int Node1_isNotConnected,
                        const amrex::IntVect Node2, const int Node2_isNotConnected,
                        const amrex::IntVect Node3, const int Node3_isNotConnected,
                        const amrex::IntVect Node4, const int Node4_isNotConnected,
                        const amrex::IntVect Node5, const int Node5_isNotConnected,
                        const amrex::IntVect Node6, const int /*Node6_isNotConnected*/,
                        const amrex::IntVect Node7, const int Node7_isNotConnected)
{

    if(covered == 1) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
       return Node6;

    } else if(covered == 2) {
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else {
          amrex::Abort("Invalid interpolation for Node 0 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node1_isNotConnected && Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node1_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node3_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else {
          amrex::Abort("Invalid interpolation for Node 0 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node1_isNotConnected && Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node1_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node3_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 0 with 4 covered cells.");
        }
    }
    amrex::Abort("Should never get here.");
    return Node0;
}

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 1                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node1_EB(const int covered,
                        const amrex::IntVect Node0, const int Node0_isNotConnected,
                        const amrex::IntVect Node1, const int /*Node1_isNotConnected*/,
                        const amrex::IntVect Node2, const int Node2_isNotConnected,
                        const amrex::IntVect Node3, const int Node3_isNotConnected,
                        const amrex::IntVect Node4, const int Node4_isNotConnected,
                        const amrex::IntVect Node5, const int Node5_isNotConnected,
                        const amrex::IntVect Node6, const int Node6_isNotConnected,
                        const amrex::IntVect Node7, const int /*Node7_isNotConnected*/)
{

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;

    } else if(covered == 2) {
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else {
            amrex::Abort("Invalid interpolation for Node 1 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else if(Node0_isNotConnected && Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node0_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node2_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else {
            amrex::Abort("Invalid interpolation for Node 1 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node0_isNotConnected && Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node0_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node2_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 1 with 4 covered cells.");
        }
    }
    amrex::Abort("Should never get here.");
    return Node1;
}

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 2                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node2_EB(const int covered,
                        const amrex::IntVect Node0, const int Node0_isNotConnected,
                        const amrex::IntVect Node1, const int Node1_isNotConnected,
                        const amrex::IntVect Node2, const int /*Node2_isNotConnected*/,
                        const amrex::IntVect Node3, const int Node3_isNotConnected,
                        const amrex::IntVect Node4, const int /*Node4_isNotConnected*/,
                        const amrex::IntVect Node5, const int Node5_isNotConnected,
                        const amrex::IntVect Node6, const int Node6_isNotConnected,
                        const amrex::IntVect Node7, const int Node7_isNotConnected)
{

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;

    } else if(covered == 2) {
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 2 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else if(Node1_isNotConnected && Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else if(Node1_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node3_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else {
            amrex::Abort("Invalid interpolation for Node 2 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node1_isNotConnected && Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node1_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node3_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 2 with 4 covered cells.");
        }
    }

    amrex::Abort("Node2: Should never get here.");
    return Node2;
}

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 3                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node3_EB(const int covered,
                        const amrex::IntVect Node0, const int Node0_isNotConnected,
                        const amrex::IntVect Node1, const int Node1_isNotConnected,
                        const amrex::IntVect Node2, const int Node2_isNotConnected,
                        const amrex::IntVect Node3, const int /*Node3_isNotConnected*/,
                        const amrex::IntVect Node4, const int Node4_isNotConnected,
                        const amrex::IntVect Node5, const int /*Node5_isNotConnected*/,
                        const amrex::IntVect Node6, const int Node6_isNotConnected,
                        const amrex::IntVect Node7, const int Node7_isNotConnected)
{

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;

    } else if(covered == 2) {
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 3 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else if(Node0_isNotConnected && Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else if(Node0_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node2_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else {
            amrex::Abort("Invalid interpolation for Node 3 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node0_isNotConnected && Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node0_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node2_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 3 with 4 covered cells.");
        }
    }

    amrex::Abort("Node3: Should never get here.");
    return Node3;
}

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 4                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node4_EB(const int covered,
                        const amrex::IntVect Node0, const int Node0_isNotConnected,
                        const amrex::IntVect Node1, const int Node1_isNotConnected,
                        const amrex::IntVect Node2, const int /*Node2_isNotConnected*/,
                        const amrex::IntVect Node3, const int Node3_isNotConnected,
                        const amrex::IntVect Node4, const int /*Node4_isNotConnected*/,
                        const amrex::IntVect Node5, const int Node5_isNotConnected,
                        const amrex::IntVect Node6, const int Node6_isNotConnected,
                        const amrex::IntVect Node7, const int Node7_isNotConnected)
{

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;

    } else if(covered == 2) {
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 4 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node0_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else if(Node0_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node5_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else {
            amrex::Abort("Invalid interpolation for Node 4 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node0_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node0_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node5_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 4 with 4 covered cells.");
        }
    }

    amrex::Abort("Node4: Should never get here.");
    return Node4;
}

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 5                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node5_EB(const int covered,
                        const amrex::IntVect Node0, const int Node0_isNotConnected,
                        const amrex::IntVect Node1, const int Node1_isNotConnected,
                        const amrex::IntVect Node2, const int Node2_isNotConnected,
                        const amrex::IntVect Node3, const int /*Node3_isNotConnected*/,
                        const amrex::IntVect Node4, const int Node4_isNotConnected,
                        const amrex::IntVect Node5, const int /*Node5_isNotConnected*/,
                        const amrex::IntVect Node6, const int Node6_isNotConnected,
                        const amrex::IntVect Node7, const int Node7_isNotConnected)
{

    if(covered == 1) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;

    } else if(covered == 2) {
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 5 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else if(Node1_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else if(Node1_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node4_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else {
            amrex::Abort("Invalid interpolation for Node 5 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node1_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node1_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node4_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 5 with 4 covered cells.");
        }
    }
    amrex::Abort("Node5: Should never get here.");
    return Node5;
}

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 6                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node6_EB(const int covered,
                        const amrex::IntVect Node0, const int /*Node0_isNotConnected*/,
                        const amrex::IntVect Node1, const int Node1_isNotConnected,
                        const amrex::IntVect Node2, const int Node2_isNotConnected,
                        const amrex::IntVect Node3, const int Node3_isNotConnected,
                        const amrex::IntVect Node4, const int Node4_isNotConnected,
                        const amrex::IntVect Node5, const int Node5_isNotConnected,
                        const amrex::IntVect Node6, const int /*Node6_isNotConnected*/,
                        const amrex::IntVect Node7, const int Node7_isNotConnected)
{

    if(covered == 1){
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;

    } else if(covered == 2) {
        if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else {
            amrex::Abort("Invalid interpolation for Node 6 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node1_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node2_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;
        } else if(Node2_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node5_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else {
            amrex::Abort("Invalid interpolation for Node 6 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node2_isNotConnected && Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 7." << std::endl;
        return Node7;
        } else if(Node2_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node5_isNotConnected && Node7_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else {
            amrex::Abort("Invalid interpolation for Node 6 with 4 covered cells.");
        }
    }

    amrex::Abort("Node6: Should never get here.");
    return Node6;
}

/*---------------------------------------------------------------------------------*
*                                                                                  *
*                                    NODE 7                                        *
*                                                                                  *
*----------------------------------------------------------------------------------*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect Node7_EB(const int covered,
                        const amrex::IntVect Node0, const int Node0_isNotConnected,
                        const amrex::IntVect Node1, const int /*Node1_isNotConnected*/,
                        const amrex::IntVect Node2, const int Node2_isNotConnected,
                        const amrex::IntVect Node3, const int Node3_isNotConnected,
                        const amrex::IntVect Node4, const int Node4_isNotConnected,
                        const amrex::IntVect Node5, const int Node5_isNotConnected,
                        const amrex::IntVect Node6, const int Node6_isNotConnected,
                        const amrex::IntVect Node7, const int /*Node7_isNotConnected*/)
{

    if(covered == 1) {
        //amrex::Print() << "Using EB in Node 1." << std::endl;
        return Node1;

    } else if(covered == 2) {
        if(Node3_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else if(Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else {
            amrex::Abort("Invalid interpolation for Node 7 with 2 covered cells.");
        }

    } else if(covered == 3){
        if(Node0_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node2_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node5_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else if(Node3_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 0." << std::endl;
        return Node0;
        } else if(Node3_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 2." << std::endl;
        return Node2;
        } else if(Node4_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 5." << std::endl;
        return Node5;
        } else {
            amrex::Abort("Invalid interpolation for Node 7 with 3 covered cells.");
        }

    } else if(covered == 4){
        if(Node3_isNotConnected && Node4_isNotConnected) {
        //amrex::Print() << "Using EB in Node 6." << std::endl;
        return Node6;
        } else if(Node3_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 4." << std::endl;
        return Node4;
        } else if(Node4_isNotConnected && Node6_isNotConnected) {
        //amrex::Print() << "Using EB in Node 3." << std::endl;
        return Node3;
        } else {
            amrex::Abort("Invalid interpolation for Node 7 with 4 covered cells.");
        }
    }

    amrex::Abort("Node7: Should never get here.");
    return Node7;
}


#endif
