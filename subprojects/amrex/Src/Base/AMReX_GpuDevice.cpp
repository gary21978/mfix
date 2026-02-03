
#include <AMReX_Arena.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Machine.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#ifdef AMREX_USE_HYPRE
#  include <_hypre_utilities.h>
#endif

#include <iostream>
#include <map>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <exception>

#if defined(AMREX_USE_CUDA)
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#if defined(AMREX_PROFILING) || defined (AMREX_TINY_PROFILING)
#if __has_include(<nvtx3/nvtx3.hpp>)
#  include <nvtx3/nvtx3.hpp>
#elif __has_include(<nvtx3/nvToolsExt.h>)
#  include <nvtx3/nvToolsExt.h>
#else
#  include <nvToolsExt.h>
#endif
#endif
#endif


#ifdef AMREX_USE_ACC
#include <openacc.h>

extern "C" {
    void amrex_initialize_acc (int);
    void amrex_finalize_acc ();
    void amrex_set_acc_stream (int);
}
#endif


namespace amrex::Gpu {

int Device::device_id = 0;
int Device::num_devices_used = 0;
int Device::num_device_partners = 1;
int Device::verbose = 0;
#ifdef AMREX_USE_GPU
int Device::max_gpu_streams = 4;
#else
int Device::max_gpu_streams = 1;
#endif

#ifdef AMREX_USE_GPU
dim3 Device::numThreadsMin      = dim3(1, 1, 1);
dim3 Device::numThreadsOverride = dim3(0, 0, 0);
dim3 Device::numBlocksOverride  = dim3(0, 0, 0);
unsigned int Device::max_blocks_per_launch = 2560;

Vector<StreamManager>   Device::gpu_stream_pool;
Vector<int>             Device::gpu_stream_index;
gpuDeviceProp_t         Device::device_prop;
int                     Device::memory_pools_supported = 0;

constexpr int Device::warp_size;


namespace {

#if defined(__CUDACC__)
    AMREX_GPU_GLOBAL void emptyKernel() {}
#endif

    void InitializeGraph(int graph_size)
    {
        amrex::ignore_unused(graph_size);

#if defined(__CUDACC__) && defined(AMREX_USE_CUDA)

        BL_PROFILE("InitGraph");

        int streams = Gpu::Device::numGpuStreams();
        cudaGraphExec_t graphExec{};
        for (int n=0; n<(graph_size); ++n)
        {
            Gpu::Device::startGraphRecording((n == 0), NULL, NULL, 0);

            // ..................
            Gpu::Device::setStreamIndex(n%streams);
            emptyKernel<<<1, 1, 0, Gpu::gpuStream()>>>();
            // ..................

            graphExec = Gpu::Device::stopGraphRecording((n == (graph_size-1)));
        }
        AMREX_CUDA_SAFE_CALL(cudaGraphExecDestroy(graphExec));
#endif
    }
}

[[nodiscard]] gpuStream_t&
StreamManager::get () {
    return m_stream;
}

void
StreamManager::sync () {
    decltype(m_free_wait_list) new_empty_wait_list{};

    {
        // lock mutex before accessing and modifying member variables
        std::lock_guard<std::mutex> lock(m_mutex);
        m_free_wait_list.swap(new_empty_wait_list);
    }
    // unlock mutex before stream sync and memory free
    // to avoid deadlocks from the CArena mutex

    // actual stream sync
    AMREX_CUDA_SAFE_CALL(cudaStreamSynchronize(m_stream));

    // synconizing the stream may have taken a long time and
    // there may be new kernels launched already, so we free memory
    // according to the state from before the stream was synced

    for (auto [arena, mem] : new_empty_wait_list) {
        arena->free(mem);
    }
}

void
StreamManager::free_async (Arena* arena, void* mem) {
    if (arena->isDeviceAccessible()) {
        std::size_t free_wait_list_size = 0;
        {
            // lock mutex before accessing and modifying member variables
            std::lock_guard<std::mutex> lock(m_mutex);
            m_free_wait_list.emplace_back(arena, mem);
            free_wait_list_size = m_free_wait_list.size();
        }
        // Limit the number of memory allocations in m_free_wait_list
        // in case the stream is never synchronized
        if (free_wait_list_size > 100) {
            sync();
        }
    } else {
        arena->free(mem);
    }
}

std::size_t
StreamManager::wait_list_size () {
    // lock mutex before accessing member variables
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_free_wait_list.size();
}

#endif

void
Device::Initialize (bool minimal, int a_device_id)
{
    amrex::ignore_unused(minimal, a_device_id);
#ifdef AMREX_USE_GPU

#if defined(AMREX_USE_CUDA) && (defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING))
    // Wrap cuda init to identify it appropriately in nvvp.
    // Note: first substantial cuda call may cause a lengthy
    // cuda API and cuda driver API initialization that will
    // be captured by the profiler. It a necessary, system
    // dependent step that is unavoidable.
    nvtxRangePush("initialize_device");
#endif

    ParmParse ppamrex("amrex");
    ppamrex.queryAdd("max_gpu_streams", max_gpu_streams);
    max_gpu_streams = std::min(max_gpu_streams, AMREX_GPU_MAX_STREAMS);
    max_gpu_streams = std::max(max_gpu_streams, 1);

    ParmParse pp("device");
    if (! pp.query("verbose", "v", verbose)) {
        pp.add("verbose", verbose);
    }

    if (amrex::Verbose()) {
        amrex::Print() << "Initializing CUDA...\n";
    }

    // Count the number of GPU devices.
    int gpu_device_count = 0;
    AMREX_CUDA_SAFE_CALL(cudaGetDeviceCount(&gpu_device_count));
    if (gpu_device_count <= 0) {
        amrex::Abort("No GPU device found");
    }

    // Now, assign ranks to GPUs. If we only have one GPU,
    // or only one MPI rank, this is easy. Otherwise, we
    // need to do a little more work.

    int n_local_procs = 1;
    amrex::ignore_unused(n_local_procs);

    if (minimal) {
        device_id = 0;
        AMREX_CUDA_SAFE_CALL(cudaGetDevice(&device_id));
    } else if (a_device_id >= 0) {
        device_id = a_device_id;
    } else if (ParallelDescriptor::NProcs() == 1) {
        device_id = 0;
    }
    else if (gpu_device_count == 1) {
        device_id = 0;
    }
    else {
        if (ParallelDescriptor::NProcsPerNode() == gpu_device_count) {
            device_id = ParallelDescriptor::MyRankInNode();
        } else if (ParallelDescriptor::NProcsPerProcessor() == gpu_device_count) {
            device_id = ParallelDescriptor::MyRankInProcessor();
        } else {
            device_id = ParallelDescriptor::MyProc() % gpu_device_count;
        }
    }

    if (gpu_device_count > 1 && ! minimal && a_device_id < 0) {
        if (Machine::name() == "nersc.perlmutter") {
            // The CPU/GPU mapping on perlmutter has the reverse order.
            device_id = gpu_device_count - device_id - 1;
            if (amrex::Verbose()) {
                amrex::Print() << "Multiple GPUs are visible to each MPI rank. Fixing GPU assignment for Perlmutter according to heuristics.\n";
            }
        } else if (Machine::name() == "olcf.frontier") {
            // The CPU/GPU mapping on fronter is documented at
            // https://docs.olcf.ornl.gov/systems/frontier_user_guide.html
            if (gpu_device_count == 8) {
                constexpr std::array<int,8> gpu_order = {4,5,2,3,6,7,0,1};
                device_id = gpu_order[device_id];
                if (amrex::Verbose()) {
                    amrex::Print() << "Multiple GPUs are visible to each MPI rank. Fixing GPU assignment for Frontier according to heuristics.\n";
                }
            }
        } else {
            if (amrex::Verbose() && ParallelDescriptor::IOProcessor()) {
                amrex::Warning("Multiple GPUs are visible to each MPI rank. This is usually not an issue. But this may lead to incorrect or suboptimal rank-to-GPU mapping.");
            }
        }
    }

    if ( ! minimal) {
        AMREX_CUDA_SAFE_CALL(cudaSetDevice(device_id));
    }

#ifdef AMREX_USE_ACC
    amrex_initialize_acc(device_id);
#endif

    initialize_gpu(minimal);

    num_devices_used = ParallelDescriptor::NProcs();

#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::NProcs() > 1 && ! minimal) {
        constexpr int len = 16;
    static_assert(std::is_same<decltype(cudaUUID_t::bytes), char[len]>());
        std::vector<char> buf(ParallelDescriptor::NProcs()*len);
        char* pbuf = buf.data();
        auto const& uuid = device_prop.uuid;
        char const* sbuf = uuid.bytes;
        MPI_Allgather(sbuf, len, MPI_CHAR, pbuf, len, MPI_CHAR,
                      ParallelDescriptor::Communicator());
        std::map<std::string,int> uuid_counts;
        std::string my_uuid;
        for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
            std::string iuuid(pbuf+i*len, len);
            if (i == ParallelDescriptor::MyProc()) {
                my_uuid = iuuid;
            }
            ++uuid_counts[iuuid];
        }
        num_devices_used = uuid_counts.size();
        num_device_partners = uuid_counts[my_uuid];

        AMREX_ALWAYS_ASSERT(num_device_partners > 0);
    }
#endif /* AMREX_USE_MPI */

        if (amrex::Verbose() && ! minimal) {
        amrex::Print() << "CUDA"
                   << " initialized with " << num_devices_used
                       << ((num_devices_used == 1) ? " device.\n"
                                                   : " devices.\n");
        if (num_devices_used < ParallelDescriptor::NProcs() && ParallelDescriptor::IOProcessor()) {
            amrex::Warning("There are more MPI processes than the number of unique GPU devices. This is not necessarily a problem.\n"
                           "For example this could happen when a device such as MI300A is partitioned into multiple subdevices.");
        }
    }

#if defined(AMREX_USE_CUDA) && (defined(AMREX_PROFILING) || defined(AMREX_TINY_PROFILING))
    nvtxRangePop();
#endif

    Device::profilerStart();

#endif /* AMREX_USE_GPU */
}

void
Device::Finalize ()
{
#ifdef AMREX_USE_GPU
    streamSynchronizeAll();
    Device::profilerStop();

    for (int i = 0; i < max_gpu_streams; ++i)
    {
        AMREX_CUDA_SAFE_CALL(cudaStreamDestroy(gpu_stream_pool[i].get()));
    }

    gpu_stream_index.clear();

#ifdef AMREX_USE_ACC
    amrex_finalize_acc();
#endif

#endif
}

void
Device::initialize_gpu (bool minimal)
{
    amrex::ignore_unused(minimal);

#ifdef AMREX_USE_GPU

    if (gpu_stream_pool.size() != max_gpu_streams) {
        // no copy/move constructor for std::mutex
        gpu_stream_pool = Vector<StreamManager>(max_gpu_streams);
    }

#ifdef AMREX_USE_CUDA
    AMREX_CUDA_SAFE_CALL(cudaGetDeviceProperties(&device_prop, device_id));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(device_prop.major >= 4 || (device_prop.major == 3 && device_prop.minor >= 5),
                                     "Compute capability must be >= 3.5");

#ifdef AMREX_GPU_STREAM_ALLOC_SUPPORT
    cudaDeviceGetAttribute(&memory_pools_supported, cudaDevAttrMemoryPoolsSupported, device_id);
#endif

#if (__CUDACC_VER_MAJOR__ < 12) || ((__CUDACC_VER_MAJOR__ == 12) && (__CUDACC_VER_MINOR__ < 4))
    if ( ! minimal ) {
        if (sizeof(Real) == 8) {
            AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
        } else if (sizeof(Real) == 4) {
            AMREX_CUDA_SAFE_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte));
        }
    }
#endif

    for (int i = 0; i < max_gpu_streams; ++i) {
        AMREX_CUDA_SAFE_CALL(cudaStreamCreate(&gpu_stream_pool[i].get()));
#ifdef AMREX_USE_ACC
        acc_set_cuda_stream(i, gpu_stream_pool[i].get());
#endif
    }

#endif

    gpu_stream_index.resize(OpenMP::get_max_threads(), 0);

    ParmParse pp("device");

    int nx = 0;
    int ny = 0;
    int nz = 0;

    pp.query("numThreads.x", nx);
    pp.query("numThreads.y", ny);
    pp.query("numThreads.z", nz);

    numThreadsOverride.x = (int) nx;
    numThreadsOverride.y = (int) ny;
    numThreadsOverride.z = (int) nz;

    nx = 0;
    ny = 0;
    nz = 0;

    pp.query("numBlocks.x", nx);
    pp.query("numBlocks.y", ny);
    pp.query("numBlocks.z", nz);

    numBlocksOverride.x = (int) nx;
    numBlocksOverride.y = (int) ny;
    numBlocksOverride.z = (int) nz;

    // Graph initialization
    int graph_init = 0;
    int graph_size = 10000;
    pp.query("graph_init", graph_init);
    pp.query("graph_init_nodes", graph_size);

    if (graph_init)
    {
        GraphSafeGuard gsg(true);
        InitializeGraph(graph_size);
    }

    max_blocks_per_launch = 4 * numMultiProcessors() * maxThreadsPerMultiProcessor() / AMREX_GPU_MAX_THREADS;

#endif
}

int
Device::deviceId () noexcept
{
    return device_id;
}

int
Device::numDevicesUsed () noexcept
{
    return num_devices_used;
}

int Device::numDevicePartners () noexcept
{
    return num_device_partners;
}

#ifdef AMREX_USE_GPU
int
Device::streamIndex (gpuStream_t s) noexcept
{
    const int N = gpu_stream_pool.size();
    for (int i = 0; i < N ; ++i) {
        if (gpu_stream_pool[i].get() == s) {
            return i;
        }
    }
    return N;
}
#endif

void
Device::setStreamIndex (int idx) noexcept
{
    amrex::ignore_unused(idx);
#ifdef AMREX_USE_GPU
    gpu_stream_index[OpenMP::get_thread_num()] = idx % max_gpu_streams;
#ifdef AMREX_USE_ACC
    amrex_set_acc_stream(idx % max_gpu_streams);
#endif
#endif
}

#ifdef AMREX_USE_GPU
gpuStream_t
Device::resetStream () noexcept
{
    gpuStream_t r = gpuStream();
    gpu_stream_index[OpenMP::get_thread_num()] = 0;
    return r;
}

gpuStream_t
Device::setStream (gpuStream_t s) noexcept
{
    gpuStream_t r = gpuStream();
    gpu_stream_index[OpenMP::get_thread_num()] = streamIndex(s);
    return r;
}
#endif

void
Device::synchronize () noexcept
{
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaDeviceSynchronize());
#endif
}

void
Device::streamSynchronize () noexcept
{
#ifdef AMREX_USE_GPU
    gpu_stream_pool[gpu_stream_index[OpenMP::get_thread_num()]].sync();
#endif
}

void
Device::streamSynchronizeAll () noexcept
{
#ifdef AMREX_USE_GPU
    for (auto& s : gpu_stream_pool) {
        s.sync();
    }
#endif
}

void
Device::freeAsync (Arena* arena, void* mem) noexcept
{
#ifdef AMREX_USE_GPU
    gpu_stream_pool[gpu_stream_index[OpenMP::get_thread_num()]].free_async(arena, mem);
#else
    arena->free(mem);
#endif
}

bool
Device::clearFreeAsyncBuffer () noexcept
{
#ifdef AMREX_USE_GPU
    bool freed_memory = false;
    for (auto& s : gpu_stream_pool) {
        if (s.wait_list_size() > 0) {
            s.sync();
            freed_memory = true;
        }
    }
    return freed_memory;
#else
    return false;
#endif
}

#if defined(__CUDACC__) && defined(AMREX_USE_CUDA)

void
Device::startGraphRecording(bool first_iter, void* h_ptr, void* d_ptr, size_t sz)
{
    if ((first_iter) && inLaunchRegion() && inGraphRegion())
    {
        // Uses passed information to do initial async memcpy in graph and
        //    links dependency to all streams using cudaEvents.

        setStreamIndex(0);
        cudaStream_t graph_stream = gpuStream();
        cudaEvent_t memcpy_event = {0};
        AMREX_CUDA_SAFE_CALL( cudaEventCreate(&memcpy_event, cudaEventDisableTiming) );

#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
        AMREX_CUDA_SAFE_CALL(cudaStreamBeginCapture(graph_stream));
#else
        AMREX_CUDA_SAFE_CALL(cudaStreamBeginCapture(graph_stream, cudaStreamCaptureModeGlobal));
#endif

        AMREX_CUDA_SAFE_CALL(cudaMemcpyAsync(d_ptr, h_ptr, sz, cudaMemcpyHostToDevice, graph_stream));
        AMREX_CUDA_SAFE_CALL(cudaEventRecord(memcpy_event, graph_stream));

        // Note: Main graph stream fixed at 0, so i starts at 1.
        //       Will need more complex logic if this changes.
        for (int i=1; i<numGpuStreams(); ++i)
        {
            setStreamIndex(i);
            AMREX_CUDA_SAFE_CALL(cudaStreamWaitEvent(gpuStream(), memcpy_event, 0));
        }
        setStreamIndex(0);

        AMREX_CUDA_SAFE_CALL( cudaEventDestroy(memcpy_event) );
    }
}

cudaGraphExec_t
Device::stopGraphRecording(bool last_iter)
{
    cudaGraphExec_t graphExec{};

    if (last_iter && inLaunchRegion() && inGraphRegion())
    {
        // Uses cudaEvents to rejoin the streams, making a single graph.
        setStreamIndex(0);
        cudaStream_t graph_stream = gpuStream();
        cudaEvent_t rejoin_event = {0};
        AMREX_CUDA_SAFE_CALL( cudaEventCreate(&rejoin_event, cudaEventDisableTiming) );

        // Note: Main graph stream fixed at 0, so i starts at 1.
        //       Will need more complex logic if this changes.
        for (int i=1; i<Gpu::Device::numGpuStreams(); ++i)
        {
            Gpu::Device::setStreamIndex(i);
            cudaEventRecord(rejoin_event, gpuStream());
            cudaStreamWaitEvent(graph_stream, rejoin_event, 0);
        }
        Gpu::Device::setStreamIndex(0);

        cudaGraph_t graph;
        AMREX_CUDA_SAFE_CALL(cudaStreamEndCapture(graph_stream, &graph));
        graphExec = instantiateGraph(graph);

        AMREX_CUDA_SAFE_CALL( cudaGraphDestroy(graph); );
        AMREX_CUDA_SAFE_CALL( cudaEventDestroy(rejoin_event) );
    }

    return graphExec;
}

cudaGraphExec_t
Device::instantiateGraph(cudaGraph_t graph)
{
    cudaGraphExec_t graphExec;

#ifdef AMREX_DEBUG
//  Implements cudaGraphInstantiate error logging feature.
//  Upon error, delays abort until message is output.
    constexpr int log_size = 1028;
    char graph_log[log_size];
    graph_log[0]='\0';

    cudaGraphInstantiate(&graphExec, graph, NULL, &(graph_log[0]), log_size);

    if (graph_log[0] != '\0')
    {
        amrex::Print() << graph_log << '\n';
        AMREX_GPU_ERROR_CHECK();
    }
#else

    AMREX_CUDA_SAFE_CALL(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));

#endif

    return graphExec;

}

void
Device::executeGraph(const cudaGraphExec_t &graphExec, bool synch)
{
    if (inLaunchRegion() && inGraphRegion())
    {
        setStreamIndex(0);
        AMREX_CUDA_SAFE_CALL(cudaGraphLaunch(graphExec, cudaStream()));
        if (synch) {
            synchronize();
        }
        resetStreamIndex();
    }
}

#endif

void
Device::mem_advise_set_preferred (void* p, std::size_t sz, int device)
{
    amrex::ignore_unused(p,sz,device);
#if defined(AMREX_USE_CUDA)
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
    {
#if defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 13)
        cudaMemLocation location = {};
        location.type = cudaMemLocationTypeDevice;
        location.id = device;
#else
        auto location = device;
#endif
        AMREX_CUDA_SAFE_CALL(
            cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, location));
    }
#endif
}

void
Device::mem_advise_set_readonly (void* p, std::size_t sz)
{
    amrex::ignore_unused(p,sz);
#if defined(AMREX_USE_CUDA)
    if (device_prop.managedMemory == 1 && device_prop.concurrentManagedAccess == 1)
    {
#if defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 13)
        cudaMemLocation location = {};
        location.type = cudaMemLocationTypeDevice;
        location.id = cudaCpuDeviceId;
#else
        auto location = cudaCpuDeviceId;
#endif
        AMREX_CUDA_SAFE_CALL(
            cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, location));
    }
#endif
}

#ifdef AMREX_USE_GPU

void
Device::setNumThreadsMin (int nx, int ny, int nz) noexcept
{
    numThreadsMin.x = nx;
    numThreadsMin.y = ny;
    numThreadsMin.z = nz;
}

void
Device::n_threads_and_blocks (const Long N, dim3& numBlocks, dim3& numThreads) noexcept
{
    numThreads = AMREX_GPU_MAX_THREADS;
    numBlocks = std::max((N + AMREX_GPU_MAX_THREADS - 1) / AMREX_GPU_MAX_THREADS, static_cast<Long>(1)); // in case N = 0
}

void
Device::c_comps_threads_and_blocks (const int* lo, const int* hi, const int comps,
                                    dim3& numBlocks, dim3& numThreads) noexcept
{
    c_threads_and_blocks(lo, hi, numBlocks, numThreads);
    numBlocks.x *= static_cast<unsigned>(comps);
}

void
Device::c_threads_and_blocks (const int* lo, const int* hi, dim3& numBlocks, dim3& numThreads) noexcept
{
    // Our threading strategy will be to allocate thread blocks
    //preferring the x direction first to guarantee coalesced accesses.
    int tile_size[] = {AMREX_D_DECL(hi[0]-lo[0]+1,hi[1]-lo[1]+1,hi[2]-lo[2]+1)};

#if (AMREX_SPACEDIM == 1)

    numThreads.x = std::min(tile_size[0], AMREX_GPU_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = 1;
    numBlocks.z = 1;

#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(tile_size[0]), AMREX_GPU_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(tile_size[1]), AMREX_GPU_MAX_THREADS / numThreads.x   );
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = (tile_size[1] + numThreads.y - 1) / numThreads.y;
    numBlocks.z = 1;

#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreads.y   ));

    numThreads.x = std::max(numThreadsMin.x, std::min(static_cast<unsigned>(tile_size[0]), numThreads.x));
    numThreads.y = std::max(numThreadsMin.y, std::min(static_cast<unsigned>(tile_size[1]), numThreads.y));
    numThreads.z = std::max(numThreadsMin.z, std::min(static_cast<unsigned>(tile_size[2]), numThreads.z));

    numBlocks.x = (tile_size[0] + numThreads.x - 1) / numThreads.x;
    numBlocks.y = (tile_size[1] + numThreads.y - 1) / numThreads.y;
    numBlocks.z = (tile_size[2] + numThreads.z - 1) / numThreads.z;

#endif

    AMREX_ASSERT(numThreads.x <= static_cast<unsigned>(device_prop.maxThreadsDim[0]));
    AMREX_ASSERT(numThreads.y <= static_cast<unsigned>(device_prop.maxThreadsDim[1]));
    AMREX_ASSERT(numThreads.z <= static_cast<unsigned>(device_prop.maxThreadsDim[2]));
    AMREX_ASSERT(numThreads.x*numThreads.y*numThreads.z <= static_cast<unsigned>(device_prop.maxThreadsPerBlock));
    AMREX_ASSERT(numThreads.x > 0);
    AMREX_ASSERT(numThreads.y > 0);
    AMREX_ASSERT(numThreads.z > 0);
    AMREX_ASSERT(numBlocks.x > 0);
    AMREX_ASSERT(numBlocks.y > 0);
    AMREX_ASSERT(numBlocks.z > 0);
}

void
Device::grid_stride_threads_and_blocks (dim3& numBlocks, dim3& numThreads) noexcept
{
    int num_SMs = device_prop.multiProcessorCount;

    int SM_mult_factor = 32;

    if (num_SMs > 0) {

        // Default to only an x loop in most cases.
        if (numThreadsMin.y == 1 && numThreadsMin.z == 1) {

            numBlocks.x = SM_mult_factor * num_SMs;
            numBlocks.y = 1;
            numBlocks.z = 1;

        } else {

            numBlocks.x = 1;
            numBlocks.y = SM_mult_factor;
            numBlocks.z = num_SMs;

        }

    } else {

        // Arbitrarily set this to a somewhat large number.

        numBlocks.x = 1000;
        numBlocks.y = 1;
        numBlocks.z = 1;

    }

#if (AMREX_SPACEDIM == 1)

    numThreads.x = std::min(device_prop.maxThreadsDim[0], AMREX_GPU_MAX_THREADS);
    numThreads.x = std::max(numThreads.x, numThreadsMin.x);
    numThreads.y = 1;
    numThreads.z = 1;

#elif (AMREX_SPACEDIM == 2)

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / numThreadsMin.y);
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / numThreads.x);
    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = 1;

#else

    numThreads.x = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[0]), AMREX_GPU_MAX_THREADS / (numThreadsMin.y * numThreadsMin.z));
    numThreads.y = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[1]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreadsMin.z));
    numThreads.z = std::min(static_cast<unsigned>(device_prop.maxThreadsDim[2]), AMREX_GPU_MAX_THREADS / (numThreads.x    * numThreads.y   ));

    numThreads.x = std::max(numThreadsMin.x, numThreads.x);
    numThreads.y = std::max(numThreadsMin.y, numThreads.y);
    numThreads.z = std::max(numThreadsMin.z, numThreads.z);

#endif

    // Allow the user to override these at runtime.

    if (numBlocksOverride.x > 0)
        numBlocks.x = numBlocksOverride.x;
    if (numBlocksOverride.y > 0)
        numBlocks.y = numBlocksOverride.y;
    if (numBlocksOverride.z > 0)
        numBlocks.z = numBlocksOverride.z;

    if (numThreadsOverride.x > 0)
        numThreads.x = numThreadsOverride.x;
    if (numThreadsOverride.y > 0)
        numThreads.y = numThreadsOverride.y;
    if (numThreadsOverride.z > 0)
        numThreads.z = numThreadsOverride.z;

}

#endif

std::size_t
Device::freeMemAvailable ()
{
#ifdef AMREX_USE_GPU
    std::size_t f;
    std::size_t t;
    AMREX_CUDA_SAFE_CALL(cudaMemGetInfo(&f,&t));
    return f;
#else
    return 0;
#endif
}

void
Device::profilerStart ()
{
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaProfilerStart());
#endif

}

void
Device::profilerStop ()
{
#ifdef AMREX_USE_CUDA
    AMREX_GPU_SAFE_CALL(cudaProfilerStop());
#endif
}

#ifdef AMREX_USE_HYPRE
void hypreSynchronize ()
{
#ifdef AMREX_USE_GPU
#if (HYPRE_RELEASE_NUMBER > 23200) || (HYPRE_RELEASE_NUMBER == 23200 && HYPRE_DEVELOP_NUMBER >= 4)
    hypre_SyncComputeStream();
#else
    hypre_SyncCudaDevice(hypre_handle()); // works for non-cuda device too
#endif
#endif
}
#endif

}
