import sys
import subprocess
import time

args = sys.argv[1:]
p = subprocess.Popen(args) #, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
print('Running: ' + ' '.join(args))
start_time = time.time()
first_time = True
while True:
    if p.poll() is not None:
        break

    ellapsed_s = time.time() - start_time
    ellapsed_min = ellapsed_s / 60
    if not first_time:
        print('Waiting for %s to finish running: %f minutes ellapsed' % 
                (args[0], ellapsed_min))
    time.sleep(5)
    first_time = False




