#pragma once

#include <mpi.h>
#include <queue>
#include <vector>

#include "include/types/checkpoint_t.h"

/**
 * Put all emerging trajectories to be calculated on a queue where every
 * process can check for new work. If all processes are finished,
 * the scheduler signals every process to finish up.
 *
 * Idea: Rank 0 solely distributes tasks and sets any process to busy.
 * Rank 0 sends its own tasks or sends to a process with available work
 * a rank id to send its own tasks to.
 */
struct task_scheduler_t{
    uint64_t current_ens;
    int current_traj;
    /**
     * A queue of all the available work to send to someone in case
     * all are busy at the moment.
     */
    std::queue<checkpoint_t> checkpoint_queue;
    /**
     * A list of all ranks with 0 being busy and 1 free for work.
     * Rank 0 is allowed to set it as busy if it sends a task or allows
     * another rank to send a task.
     * Every other rank can use MPI_Put for its own idx to signal it is free.
     * If free_worker[0] == 2, then all processes are finished.
     */
    std::vector<std::int8_t> free_worker;
    /**
     * Number of work available at rank idx but not sent yet.
     * Rank 0 can reduce it with each signal by one. Each rank can modify
     * its own idx with MPI_Put in case it just took over its own available
     * work or found new one that could not be sent.
     */
    std::vector<std::int8_t> work_available;

    /**
     * If 1 at idx then send work to process idx.
     */
    // std::vector<std::int8_t> send_list;

    /**
     * Store the maximum ensemble id that had been created.
     * rank 0: Always has the maximum number at hand.
     * rank > 0: not necessarily up to date.
     */
    uint64_t max_ensemble_id;

    MPI_Win free_window;
    MPI_Win work_window;
    MPI_Win ens_window;
    int my_rank;

    task_scheduler_t(const int &rank, const int &n_processes,
        const int &simulation_mode);

    /**
     * Check if someone is available to send a task from queue to and send it.
     * @param checkpoint on out: checkpoint if receiver and sender are
     * the same
     * @param send_to_self If true: allow to store checkpoint in
     * checkpoint for further use.
     * @return True if own queue is where it shall be sent to.
     */
    bool send_task(checkpoint_t &checkpoint, const bool send_to_self = true);

    /**
     * Checks if any worker is free and send it a new task.
     * If all workers are busy, add the task to its queue for
     * later sends.
     */
    void send_new_task(checkpoint_t &checkpoint);

    /**
     * Busy waiting until a new task is available and returns true.
     * If all processes are waiting, then no more work is available
     * and it returns false, signaling that all processes are finished.
     * Returns immediately if static scheduling is enabled.
     *
     * @return False if all processes are finished and there is no new work
     * True if new work had been found and loaded.
     */
    bool receive_task(checkpoint_t &checkpoint);

    /**
     * Broadcast to other processes that this process is free for work.
     */
    void signal_free();

    /**
     * Signal to rank 0 how much work is available.
     */
    void signal_work_avail();

    void set_n_ensembles(const int &n) {n_ensembles = n;}
    void set_n_trajectories(const int &n) {n_trajectories = n;}

 private:
    bool static_scheduling;
    uint64_t n_ensembles;
    uint64_t n_trajectories;
    uint64_t n_processes;
    /**
     * Only for rank 0. Checks for any available work on other
     * processes and signals them to send it.
     * Is called by send_task().
     */
    void signal_send_task();
};
