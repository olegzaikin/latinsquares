import sys

fname = sys.argv[1]
total_models_num = 0
models_str_num = 0
with open(fname, 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        if not line.startswith('c Models'):
            continue
        models_str_num += 1
        words = line.split()
        assert(len(words) == 4)
        total_models_num += int(words[3])
print('models in total : ' + str(total_models_num))
print('strings with models numbers in total : ' + str(models_str_num))
