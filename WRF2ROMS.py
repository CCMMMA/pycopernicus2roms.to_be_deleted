import argparse
import time

from jncregridder.roms.ROMSGrid import ROMSGrid
from jncregridder.roms.ROMSWind import ROMSWind

import warnings
warnings.filterwarnings("ignore", category=UserWarning)


class WRF2ROMS:
    def __init__(self, romsGridFilename, wrfFilename, forcingFilename):
        self.romsGridFilename = romsGridFilename
        self.wrfFilenameOrFilesPath = wrfFilename
        self.romsForcingFilename = forcingFilename

        # Open the ROMS Grid
        romsGrid = ROMSGrid(self.romsGridFilename)

        # Open ROMS wind File
        romsWind = ROMSWind(self.romsForcingFilename, romsGrid, self.wrfFilenameOrFilesPath)
        start = time.time()
        romsWind.make()
        end = time.time()
        print(end-start)

        romsWind.close()


def parser():
    parser = argparse.ArgumentParser(description="WRF2ROMS")
    parser.add_argument("--romsGridFilename", type=str, required=True, help="Path to grid copernicus-data file")
    parser.add_argument("--wrfFilename", type=str, required=True, help="Path to wrf file/files")
    parser.add_argument("--forcingFilename", type=str, required=True, help="")
    return parser


def main():
    arg_parser = parser()
    args = arg_parser.parse_args()

    WRF2ROMS(args.romsGridFilename, args.wrfFilename, args.forcingFilename)


if __name__ == '__main__':
    main()
