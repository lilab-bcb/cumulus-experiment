#!/bin/bash

function get_mem_info() {
        # /proc/meminfo
        cat /proc/meminfo
}

function get_mem_available() {
        # mem unused from /proc/meminfo
        get_mem_info | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print $2 }'
}

function get_mem_total() {
        # mem total from /proc/meminfo
        get_mem_info | grep MemTotal | awk 'BEGIN { FS=" " } ; { print $2 }'
}

function get_mem_usage() {
        # memTotal and memAvailable
        local -r mem_total=$(get_mem_total)
        local -r mem_available=$(get_mem_available)

        # usage = 100 * mem_used / mem_total
        local -r mem_used=$(($mem_total-$mem_available))
        echo "$mem_used" "$mem_total" "%"| awk '{ print 100*($1/$2)$3 }'
}

function print_usage() {
        echo [$(date)]
        echo \* Memory usage: "$(get_mem_usage)"
}

function print_summary() {
        # display header information
        echo ==================================
        echo =========== MONITORING ===========
        echo ==================================

        # summary info
        echo --- General Information ---
        # number of cores
        echo \#CPU: $(nproc)
        # multiply by 10^-6 to convert KB to GB
        echo Total Memory: $(echo $(get_mem_total) 1000000 | awk '{ print $1/$2 }')G
        echo Total Disk space: $(df -h | grep cromwell_root | awk '{ print $2}')
}

function main() {
        # disk, mem and cpu general statisitcs
        print_summary

        # sleep b/w getting usage and intially storing the cpu_previous usage values
        # this is b/c cpu usage values are time dependent
        # to calculate cpu usage, values must be determined from 2 diff time stamps
       	if [ -z "$MONITOR_SCRIPT_SLEEP" ]; then
			MONITOR_SCRIPT_SLEEP=1
		fi
        # get usage of disk, cpu and mem every MONITOR_SCRIPT_SLEEP sec
        echo
        echo --- Runtime Information ---

        sleep "$MONITOR_SCRIPT_SLEEP";
        while true; do print_usage; sleep "$MONITOR_SCRIPT_SLEEP"; done
}

main
