import re, sys

def get_max_mem(filename, memory_limit):
    max_mem_percent = 0.0
    with open(filename, 'r') as f:
        for line in f:
            if re.match('[*] Memory usage:.*', line):
                n = float(line.rstrip().split(' ')[-1][:-1])
                max_mem_percent = n if n > max_mem_percent else max_mem_percent

    return memory_limit * max_mem_percent / 100, max_mem_percent

if __name__ == '__main__':
    filename = sys.argv[1]
    max_mem, max_percent = get_max_mem(filename, memory_limit = 256)
    print("Maximum memory usage is {size:.4f} GB ({percent}%).".format(size = max_mem, percent = max_percent))