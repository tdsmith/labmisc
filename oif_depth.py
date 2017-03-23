import codecs
import ConfigParser as cp
import sys

def read_oif(fp):
    # returns a tuple (start, end, slices, step)
    parser = cp.SafeConfigParser()
    parser.readfp(fp)
    scan_mode = parser.get("Acquisition Parameters Common", "ScanMode")
    section = "Axis 3 Parameters Common"
    start = parser.getfloat(section, "StartPosition")
    end = parser.getfloat(section, "EndPosition")
    if scan_mode == '"XYZ"':
        slices = parser.getint(section, "MaxSize")
    else:
        slices = 0
    step = parser.getfloat(section, "Interval")
    units = parser.get(section, "UnitName")
    divisor = 1000. if units == '"nm"' else 1.
    return (start/divisor, end/divisor, slices, step/divisor)


def main():
    print "start\tend\tslices\tstep\tfilename"
    for i in sys.argv[1:]:
        with codecs.open(i, "r", encoding="utf-16") as f:
            start, end, slices, step = read_oif(f)
        print "%.2f\t%.2f\t%d\t%.2f\t%s" % (
            start, end, slices, step, i
        )

if __name__ == "__main__":
    main()
