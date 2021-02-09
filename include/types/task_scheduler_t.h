#pragma once

#include "include/types/checkpoint_t.h"

/**
 * Put all emerging trajectories to be calculated on a queue where every
 * process can check for new work. If all processes are finished,
 * the scheduler signals every process to finish up.
 */
struct task_scheduler_t{


    task_scheduler_t(const int &rank, const int &n_processes);
    /**
     * Busy waiting for new work. If all processes signal that they are
     * waiting, there is no more work and we return False. If new work
     * is found, load it and return True.
     *
     * @return False if all processes are finished and there is no new work
     * True if new work had been found and loaded.
     */
    bool check_queue();

    /**
     * Add a simulation to the queue using a checkpoint.
     */
    void push_back(checkpoint_t &checkpoint);

    /**
     *
     */
    void pop_front();

};