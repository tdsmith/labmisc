import numpy as np
import scipy as sp
import scipy.interpolate
from serial import Serial
from time import sleep, time
import sys

_state = []
_temp = None

mean = lambda x: float(sum(x))/len(x)

def read(s):
    global _state, _temp
    buf = s.read(512)
    while buf:
        _state.append(buf)
        buf = s.read(512)
    lines = ''.join(_state).split('\r')
    if not lines[0].startswith('sc'):
        lines = lines[1:]
    if not lines:
        return
    last = lines.pop()
    if last:
        _state = [last]
    else:
        _state = []
    lines = [[int(Vcc), int(Vr), int(switch)] for (sc, Vcc, Vr, switch) in [line.split() for line in lines]]
    mean_Vcc = mean([line[0] for line in lines]) * 10.0/2048
    mean_Vr = mean([line[1] for line in lines]) * 10.0/2048
    R = (mean_Vcc/mean_Vr - 1) * 10000.
    T = _temp(R)
    print time(), len(lines), mean_Vcc, mean_Vr, R, T
    sys.stdout.flush()

def main(tty):
    global _temp
    lookup = np.flipud(np.loadtxt('therm_lookup.csv', delimiter=','))
    _temp = sp.interpolate.interp1d(lookup[:,1] * 1000, lookup[:,0], kind='cubic')
    s = Serial(tty, baudrate=115200, rtscts=True, timeout=0)
    s.flushInput()
    s.write('\rstop\rinfo 0\r')
    s.write('\r'.join(['', 'asc', 'slist 0 x0', 'slist 1 x1', 'slist 2 x8', '']))
    s.write('start\r')
    sleep(1)
    s.flushInput()
    try:
        while True:
            sleep(1)
            read(s)
    finally:
        s.write('stop\r')
        s.close()

if __name__ == '__main__':
    import sys
    main(sys.argv[1])
