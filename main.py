import pybedtools
from pybedtools import BedTool
import argparse

def main(path1, path2, savebedpath, savecsvpath, tmpfilepath, RegInterval, Thread, Length = 4000):
    # path1: path for GEMs (i.e. ___ALL.region.PEanno)
    # path2: path for Region (i.e. ____PETcnt_G9.motifannot)
    # savebedpath: path for saving extracted GEMs in .bed
    # savecsvpath: path for saving summary table in .csv
    # tmpfilepath: path for saving tmpfiles produced by pybedtool, a directory
    # Thread: for naming the csv file. (i.e. '0')
    # Length: Length of extension. Default = 4000 (int)

    # path1: path for GEMs (i.e. ___ALL.region.PEanno)
    # path2: path for Region (i.e. ____PETcnt_G9.motifannot)
    # savebedpath: path for saving extracted GEMs in .bed
    # savecsvpath: path for saving summary table in .csv
    # tmpfilepath: path for saving tmpfiles produced by pybedtool, a directory
    # Thread: for naming the csv file. (i.e. '0')
    # Length: Length of extension. Default = 4000 (int)
    pybedtools.helpers.cleanup()
    pybedtools.set_tempdir(tmpfilepath)

    # Specify for the path of ___ALL.region.PEanno and import it (GEMs)
    # path1 = 'Minji_data/SHG0180-181-182NR_hg38_cohesin_FDR_0.1_ALL_motifext4kbboth.region.PEanno'
    ChIA_Drop = BedTool(path1)

    # Specify for the path of ____PETcnt_G9.motifannot and import it (anchors, regions)
    # path2 = 'Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot.sorted.domains'
    Region = BedTool(path2)

    # Remove unnecessary entries
    Region_short = Region.groupby(g=[1,2,6,12,14,20,8,9,16,21], c=[12], o=['count'])
    # Region_short.moveto('Region_short.bed')
    # Region_short = BedTool('Region_short.bed')
    Max_iter = Region_short.count()
    # Length = 4000


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path1',type = str)
    parser.add_argument('--path2',type = str)
    parser.add_argument('--savebedpath',type = str)
    parser.add_argument('--savecsvpath',type = str)
    parser.add_argument('--tmpfilepath',type = str)

    parser.add_argument('--Start_pos',type = int)
    parser.add_argument('--End_pos',type = int)
    parser.add_argument('--Length',type = int)
    parser.add_argument('--Thread',type = str)
    args = parser.parse_args()

    path1 = args.path1
    path2 = args.path2
    savebedpath = args.savebedpath
    savecsvpath = args.savecsvpath
    tmpfilepath = args.tmpfilepath
    Start_pos = args.Start_pos
    End_pos = args.End_pos
    Thread = args.Thread
    Length = args.Length

    RegInterval = range(Start_pos, End_pos)
    # Note that it is 0-based.
    # It will process through exactly Start_pos to End_pos-1 (i.e. range(Start_pos,End_pos))
    # path1: path for GEMs (i.e. ___ALL.region.PEanno)
    # path2: path for Region (i.e. ____PETcnt_G9.motifannot)
    # savebedpath: path for saving extracted GEMs in .bed
    # savecsvpath: path for saving summary table in .csv
    # tmpfilepath: path for saving tmpfiles produced by pybedtool, a directory
    # Thread: for naming the csv file. (i.e. '0')
    # Length: Length of extension. Default = 4000 (int)
    main(path1, path2, savebedpath, savecsvpath, tmpfilepath, RegInterval, Thread, Length)