#include "include/types/task_scheduler_t.h"


task_scheduler_t::task_scheduler_t(
    const int &rank,
    const int &n_processes,
    const int &simulation_mode) {
    if ((simulation_mode == grid_sensitivity) || (simulation_mode == trajectory_sensitivity)
        || (simulation_mode == create_train_set)) {
        static_scheduling = true;
    } else {
        static_scheduling = false;
    }
    signal_sent = false;

#ifdef COMPRESS_OUTPUT
    // We need a way to know if any process is still running for
    // collective write (and compression) of data.
    free_worker.resize(n_processes);
    if (rank == 0) {
        std::fill(free_worker.begin()+1, free_worker.end(), 1);
        free_worker[0] = 0;
    } else {
        std::fill(free_worker.begin(), free_worker.end(), 0);
    }
    work_available.resize(n_processes);
    std::fill(work_available.begin(), work_available.end(), 0);
    max_ensemble_id = 0;

    SUCCESS_OR_DIE(
        MPI_Win_create(
            free_worker.data(),
            n_processes*sizeof(int),
            sizeof(int),
            MPI_INFO_NULL,
            MPI_COMM_WORLD,
            &free_window));
    SUCCESS_OR_DIE(
        MPI_Win_create(
            work_available.data(),
            n_processes*sizeof(int),
            sizeof(int),
            MPI_INFO_NULL,
            MPI_COMM_WORLD,
            &work_window));
    SUCCESS_OR_DIE(
        MPI_Win_create(
            &max_ensemble_id,
            sizeof(std::uint64_t),
            sizeof(std::uint64_t),
            MPI_INFO_NULL,
            MPI_COMM_WORLD,
            &ens_window));
    MPI_Win_fence(0, free_window);
    MPI_Win_fence(0, work_window);
    MPI_Win_fence(0, ens_window);
    if (static_scheduling) {
        current_traj = rank-n_processes;
        current_ens = 0;
    } else {
        current_traj = 0;
        current_ens = 0;
    }
#else
    if (!static_scheduling) {
        free_worker.resize(n_processes);
        if (rank == 0) {
            std::fill(free_worker.begin()+1, free_worker.end(), 1);
            free_worker[0] = 0;
        } else {
            std::fill(free_worker.begin(), free_worker.end(), 0);
        }
        work_available.resize(n_processes);
        std::fill(work_available.begin(), work_available.end(), 0);
        max_ensemble_id = 0;

        SUCCESS_OR_DIE(
            MPI_Win_create(
                free_worker.data(),
                n_processes*sizeof(int),
                sizeof(int),
                MPI_INFO_NULL,
                MPI_COMM_WORLD,
                &free_window));
        SUCCESS_OR_DIE(
            MPI_Win_create(
                work_available.data(),
                n_processes*sizeof(int),
                sizeof(int),
                MPI_INFO_NULL,
                MPI_COMM_WORLD,
                &work_window));
        SUCCESS_OR_DIE(
            MPI_Win_create(
                &max_ensemble_id,
                sizeof(std::uint64_t),
                sizeof(std::uint64_t),
                MPI_INFO_NULL,
                MPI_COMM_WORLD,
                &ens_window));
        MPI_Win_fence(0, free_window);
        MPI_Win_fence(0, work_window);
        MPI_Win_fence(0, ens_window);
#ifdef DEBUG_WINDOW
        std::cout << rank << "created windows \n";
#endif
        current_traj = 0;
        current_ens = 0;
    } else {
        current_traj = rank-n_processes;
        current_ens = 0;
    }
#endif
    my_rank = rank;
    this->n_processes = n_processes;
}


bool task_scheduler_t::send_task(
    checkpoint_t &checkpoint,
    const bool send_to_self) {
    if (static_scheduling) return true;

    bool orig_is_dest = false;
    bool locked = false;
    if (my_rank == 0 || !checkpoint_queue.empty()) {
        SUCCESS_OR_DIE(
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, my_rank, 0, work_window));
        locked = true;
    }
    // Send a task from own queue if some worker is free or
    // if send_list already has a worker in mind
    if (!checkpoint_queue.empty()) {
        auto it = std::find(free_worker.begin(), free_worker.end(), 1);
        while (it != free_worker.end() && !checkpoint_queue.empty()) {
            auto idx = std::distance(free_worker.begin(), it);
            checkpoint = checkpoint_queue.front();
            checkpoint_queue.pop();
            if (my_rank != idx)
                checkpoint.send_checkpoint(idx);
            else
                orig_is_dest = true;

            if (my_rank == 0) {
                free_worker[idx] = 0;
                work_available[0] = checkpoint_queue.size();
            }
            it = std::find(free_worker.begin(), free_worker.end(), 1);
        }
        // in case a checkpoint is still left over and had not been sent
        // to oneself.
        if (!checkpoint_queue.empty() && !orig_is_dest && send_to_self) {
            checkpoint = checkpoint_queue.front();
            checkpoint_queue.pop();
            orig_is_dest = true;
            if (my_rank == 0) {
                // free_worker[0] = 0;
                int busy = 0;
                int current_value = 1;
                int compare_value = 1;
                SUCCESS_OR_DIE(
                    MPI_Compare_and_swap(
                        &busy,              // source
                        &compare_value,     // compare
                        &current_value,     // result
                        MPI_INT,            // datatype
                        0,                  // target rank
                        0,                  // displacement in window
                        free_window));
                work_available[0] = checkpoint_queue.size();
            }
        }
    }
    if (my_rank == 0 && !orig_is_dest) signal_send_task();
    if (locked) {
        SUCCESS_OR_DIE(MPI_Win_unlock(my_rank, work_window));
    }
    return orig_is_dest;
}


void task_scheduler_t::send_new_task(
    checkpoint_t &checkpoint) {
    if (static_scheduling) return;
    SUCCESS_OR_DIE(
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, my_rank, 0, work_window));
    auto it = std::find(free_worker.begin(), free_worker.end(), 1);
    if (it != free_worker.end()) {
        auto idx = std::distance(free_worker.begin(), it);
        checkpoint.send_checkpoint(idx);
        int busy = 0;
        int current_value = 1;
        int compare_value = 1;
        SUCCESS_OR_DIE(
            MPI_Compare_and_swap(
                &busy,              // source
                &compare_value,     // compare
                &current_value,     // result
                MPI_INT,            // datatype
                0,                  // target rank
                idx,                // displacement in window
                free_window));
    } else {
        checkpoint_queue.push(checkpoint);
        work_available[my_rank]++;
        if (my_rank != 0) {
            signal_work_avail();
        }
    }
    SUCCESS_OR_DIE(MPI_Win_unlock(my_rank, work_window));
}


void task_scheduler_t::signal_send_task() {
    if (static_scheduling) return;
    // Check for any remaining tasks on any other process and signal
    // them where to send it.
    if (my_rank == 0) {
        uint32_t work_idx = 0;
        for (auto &worker : free_worker) {
            if (worker == 1) {
                for (uint32_t i=work_idx; i < work_available.size(); ++i) {
                    if (work_available[i] > 0) {
                        // signal to work_idx to send to worker
                        int free = 1;
                        int current_value = 0;
                        int compare_value = 0;
                        SUCCESS_OR_DIE(
                            MPI_Compare_and_swap(
                                &free,              // source
                                &compare_value,     // compare
                                &current_value,     // result
                                MPI_INT,            // datatype
                                work_idx,           // target rank
                                worker,             // displacement in window
                                free_window));

                        // Mark worker as busy and remove work
                        int busy = 0;
                        current_value = 1;
                        compare_value = 1;
                        SUCCESS_OR_DIE(
                            MPI_Compare_and_swap(
                                &busy,              // source
                                &compare_value,     // compare
                                &current_value,     // result
                                MPI_INT,            // datatype
                                my_rank,            // target rank
                                worker,             // displacement in window
                                free_window));
                        work_available[i]--;
                        work_idx = i;
                        break;  // continue with next free worker
                    }
                }
            }
        }
    }
}


bool task_scheduler_t::all_free() {
    // Only rank 0 can determine if every process is finished.
    if (my_rank == 0) {
        auto it = std::find(free_worker.begin(), free_worker.end(), 0);
        if (it == free_worker.end()) {
             // signal everyone that we are finished somehow
            free_worker[0] = 2;
            // We send it to every worker, however broadcast a single
            // value everytime we are here is not a good idea.
            for (uint32_t i=1; i < free_worker.size(); ++i) {
                int current_value = 0;
                int compare_value = 0;
                SUCCESS_OR_DIE(
                    MPI_Compare_and_swap(
                        free_worker.data(),  // source
                        &compare_value,     // compare
                        &current_value,     // result
                        MPI_INT,         // datatype
                        i,                  // target rank
                        0,            // displacement in window
                        free_window));
            }
            return true;
        }
    } else {
        if (free_worker[0] == 2) {
            return true;
        }
    }
    return false;
}

bool task_scheduler_t::receive_task(
    checkpoint_t &checkpoint) {
    if (static_scheduling) {
        uint64_t next_idx = current_traj + n_processes;
        if (next_idx >= n_trajectories) {
            next_idx = next_idx%n_trajectories;
            current_ens++;
        }
#ifdef COMPRESS_OUTPUT
        if (current_ens >= n_ensembles) {
            // We are finished and can signal that to all processes
            signal_free();
            return false;
        }
#else
        if (current_ens >= n_ensembles) return false;
#endif
        current_traj = next_idx;
        return true;
    }
    // check if work is directly available
    int all_done = 0;
    bool task_sent = false;
    while (!checkpoint.receive_checkpoint()) {
        // check if work is available on own queue ?
        // Send any tasks in case it didn't happen yet
        if (send_task(checkpoint)) task_sent = true;  // return true;
        // if not, signal free to the master and continue busy waiting
        if (!signal_sent && my_rank != 0 && !task_sent) {
            signal_free();
            signal_sent = true;
        }
        // if all are free, return
        if (my_rank == 0 && !task_sent) {
            SUCCESS_OR_DIE(
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, my_rank, 0, work_window));
            int free = 1;
            int current_value = 0;
            int compare_value = 0;
            SUCCESS_OR_DIE(
                MPI_Compare_and_swap(
                    &free,              // source
                    &compare_value,     // compare
                    &current_value,     // result
                    MPI_INT,            // datatype
                    0,                  // target rank
                    my_rank,            // displacement in window
                    free_window));

            auto it = std::find(free_worker.begin(), free_worker.end(), 0);
            bool no_more_work = true;
            for (const auto &work : work_available) {
                if (work > 0) {
                    no_more_work = false;
                    break;
                }
            }
            // Get that work done by someone
            if (!no_more_work)
                signal_send_task();
            if (it == free_worker.end() && no_more_work) {
                all_done = 1;
            }
            SUCCESS_OR_DIE(MPI_Win_unlock(my_rank, work_window));
        }
        SUCCESS_OR_DIE(MPI_Bcast(&all_done, 1, MPI_INT, 0, MPI_COMM_WORLD));
        if (task_sent) return true;
        if (all_done == 1) return false;
    }
    signal_sent = false;

    if (all_done == 1) return false;
    return true;
}


void task_scheduler_t::signal_free() {
    // Set to free (1) if current status is busy (0).
    // One attempt is enough. Otherwise we could see an overlap
    // between sending data to this process and signaling free in
    // rare cases.
    // Or in other words: If process 0 already thinks this process is
    // free (1) then there is no need for a swap.
    int free = 1;
    int current_value = 0;
    int compare_value = 0;
    SUCCESS_OR_DIE(
        MPI_Compare_and_swap(
            &free,              // source
            &compare_value,     // compare
            &current_value,     // result
            MPI_INT,            // datatype
            0,                  // target rank
            my_rank,            // displacement in window
            free_window));
}


void task_scheduler_t::signal_work_avail() {
    SUCCESS_OR_DIE(
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, work_window));
    int work = checkpoint_queue.size();
    SUCCESS_OR_DIE(
        MPI_Put(
            &work,          // source
            1,              // count
            MPI_INT,        // datatype
            0,              // target
            my_rank,        // displacement in target
            1,              // target count
            MPI_INT,        // target datatype
            work_window));
    SUCCESS_OR_DIE(MPI_Win_unlock(0, work_window));
}
