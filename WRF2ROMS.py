import argparse
import os

from jncregridder.roms.ROMSGrid import ROMSGrid
from jncregridder.roms.ROMSWind import ROMSWind
from jncregridder.wrf.WRFData import WRFData

import warnings
warnings.filterwarnings("ignore", category=UserWarning)


class WRF2ROMS:
    def __init__(self, romsGridFilename, wrfFilename, timeOffset, forcingFilename):
        self.romsGridFilename = romsGridFilename
        self.wrfFilenameOrFilesPath = wrfFilename
        self.timeOffset = timeOffset
        self.romsForcingFilename = forcingFilename

        # Open the ROMS Grid
        romsGrid = ROMSGrid(self.romsGridFilename)

        # Open ROMS wind copernicus-data
        romsWind = ROMSWind(self.romsForcingFilename, romsGrid)

        folder = os.path.abspath(self.wrfFilenameOrFilesPath)
        if os.path.isdir(folder):
            listOfFiles = sorted(os.listdir(folder))
            for filename in listOfFiles:
                filepath = os.path.join(folder, filename)
                if os.path.isfile(filepath) and filename.startswith("wrf"):
                    print("File", filepath)
                    wrfData = WRFData(filepath)
                    romsWind.add(wrfData, timeOffset)
        else:
            print(self.wrfFilenameOrFilesPath)
            wrfData = WRFData(self.wrfFilenameOrFilesPath)
            romsWind.add(wrfData, timeOffset)

        romsWind.close()


def parser():
    parser = argparse.ArgumentParser(description="WRF2ROMS")
    parser.add_argument("--romsGridFilename", type=str, required=True, help="Path to grid copernicus-data file")
    parser.add_argument("--wrfFilename", type=str, required=True, help="Path to wrf file/files")
    parser.add_argument("--timeOffset", type=int, required=True, help="")
    parser.add_argument("--forcingFilename", type=str, required=True, help="")
    return parser


def main():
    arg_parser = parser()
    args = arg_parser.parse_args()

    WRF2ROMS(args.romsGridFilename, args.wrfFilename, args.timeOffset, args.forcingFilename)


if __name__ == '__main__':
    main()
