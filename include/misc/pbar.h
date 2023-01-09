#pragma once

#include <math.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include <chrono>  // NOLINT [build/c++11]
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class ProgressBar {
 public:
    ProgressBar() {
        end_step = 0;
        update_every = 1;
        description = "simulation step";
        unit = "Hz";
        out = nullptr;
        desc_width = description.length();
        current_step = 0;
    }
    ProgressBar(
        const uint64_t end_step_,
        const uint64_t update_every_,
        const std::string description_,
        std::ostream& out_) {
        end_step = end_step_;
        end_step_string = std::to_string(end_step);
        update_every = (update_every_ > end_step) ? end_step : update_every_;
        description = description_;
        unit = "Hz";
        out = &out_;
        // Description width with spaces, additional characters and some additional space
        desc_width = description.length() + unit.length()*2 + 16 + 40 + 14;
        t_first = std::chrono::system_clock::now();
        current_step = 0;
    }

    void reset() {
        t_first = std::chrono::system_clock::now();
        desc_width = description.length() + unit.length()*2 + 16 + 40 + 14;
    }
    void finish() {
        if (update_every == 0) return;
        current_step = end_step-1;
        progress();
        *out << "\n" << std::flush;
    }
    void set_steps(const uint64_t end_step_,
        const uint64_t update_every_) {
        end_step = end_step_;
        end_step_string = std::to_string(end_step);
        update_every = (update_every_ > end_step) ? end_step : update_every_;
    }
    void set_outstream(
        std::ostream& out_) {
        out = &out_;
    }
    void set_current_step(const uint64_t current_step_) {
        current_step = current_step_;
    }
    void progress() {
        // While debugging, the bar is not useful.
#if !defined(TRACE_SAT) && !defined(TRACE_ENV) && !defined(TRACE_QV) && !defined(TRACE_QC) && !defined(TRACE_QR)
#if !defined(TRACE_QS) && !defined(TRACE_QI) && !defined(TRACE_QG) && !defined(TRACE_QH)
        current_step++;
        if (update_every == 0) return;
        if (current_step%update_every != 0 && current_step != end_step) return;
        if (current_step == 0) return;
        if (current_step > end_step)
            end_step = current_step + 1;
        auto now = std::chrono::system_clock::now();
        double dt_total = ((std::chrono::duration<double>)(now - t_first)).count();
        // Get total time string
        std::string time = " in ";
        if (dt_total >= 3600)
            time = time + std::to_string(static_cast<int>(dt_total/3600)) + "h ";
        if (dt_total >= 60)
            time = time + std::to_string(static_cast<int>(dt_total)%3600/60) + "min ";

        time = time + std::to_string(static_cast<int>(dt_total)%60) + "s ";
        int window_width = get_console_width();
        // Get average amount of steps per second
        double dt_step = round(current_step/dt_total);
        std::stringstream right_side;
        right_side << std::fixed << std::setprecision(3) << description
                   << ": " << current_step << "/"
                   << end_step_string << time << "(";
        if (dt_step > 1e6) {
            right_side << static_cast<int>(dt_step/1e6) << "MHz)";
        } else if (dt_step > 1e3) {
            right_side << static_cast<int>(dt_step/1e3) << "kHz)";
        } else {
            right_side << static_cast<int>(dt_step) << "Hz)";
        }
        std::string right_string = right_side.str();

        // Get estimated remaining time "Rem. xxmin xxs"
        uint64_t rem_time = (end_step - current_step)/dt_step;
        std::string rem_string = " Rem. ";
        if (rem_time >= 3600)
            rem_string = rem_string + std::to_string(static_cast<int>(rem_time/3600)) + "h ";
        if (rem_time >= 60)
            rem_string = rem_string + std::to_string(static_cast<int>(rem_time%3600/60)) + "min ";
        rem_string = rem_string + std::to_string(static_cast<int>(rem_time%60)) + "s";

        // Check how long the bar can be
        if (right_string.length()+1+rem_string.length() > desc_width) {
            desc_width = right_string.length() + rem_string.length();
        } else {
            while (right_string.length()+rem_string.length() < desc_width) {
                rem_string = rem_string + " ";
            }
        }
        // Get the progressbar
        int bar_width_max = window_width-desc_width-1;
        double current_box = static_cast<double>(current_step)/static_cast<double>(end_step) * bar_width_max;
        int n_full = current_box;
        std::string bar;
        for (int i=0; i < n_full; i++) bar += bars[8];
        // fill the one inbetween
        int bar_width = n_full;
        if (current_box > n_full) {
            bar += bars[static_cast<uint32_t>(bars.size()*(current_box-n_full))];
            bar_width++;
        }
        // the rest
        for (int i=0; i < bar_width_max-bar_width; i++) bar += bars[0];
        bar += right_pad;

        *out << bar << right_string << rem_string << "\r" << std::flush;
#endif
#endif
    }

 private:
    std::vector<std::string> bars = {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"};
    std::chrono::time_point<std::chrono::system_clock> t_first = std::chrono::system_clock::now();
    uint64_t end_step;
    std::string end_step_string;
    uint64_t update_every;
    const std::string right_pad = "▏";
    std::string description;
    std::string unit;
    uint64_t desc_width;
    std::ostream* out;
    uint64_t current_step;

    int get_console_width() {
        struct winsize win;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &win);
        // Weird numbers could mean, this is running using a windows
        // terminal, ie. via ssh.
        if (win.ws_col > 5000 || win.ws_col < 150) return 150;
        return win.ws_col;
    }
};
