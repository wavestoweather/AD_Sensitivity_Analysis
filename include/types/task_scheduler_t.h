#pragma once

#include <mpi.h>
#include <queue>

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

    MPI_Win free_window;
    MPI_Win work_window;
    int my_rank;



    task_scheduler_t(const int &rank, const int &n_processes);

    /**
     * Check if someone is available to send a task to and send it.
     */
    void send_task();

    /**
     * Checks if any worker is free and send it one task.
     * If all workers are busy, add the task to its queue for
     * later sends.
     */
    void send_task(checkpoint_t &checkpoint);

    /**
     * Busy waiting until a new task is available and returns true.
     * If all processes are waiting, then no more work is available
     * and it returns false, signaling that all processes are finished.
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
  private:
    /**
     * Only for rank 0. Checks for any available work on other
     * processes and signals them to send it.
     * Is called by send_task().
     */
    void signal_send_task();
};