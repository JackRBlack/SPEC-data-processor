########## INFO ##########
print("########################################")
print("Project: spec_data_analyser")
print("Version: 0.8.3")
print("Last Update: 2019.08.28")
print("----------------------------------------")
print("Author: Wenjie Chen")
print("E-mail: wenjiechen@pku.edu.cn")
print("########################################")
##########################

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

import csv
import time

import os
import shutil

from IPython.display import Video

from UB_matrix import *

class specdata:
    def __init__(self, NAME="NULL"):
        self.PROJECT_NAME = NAME

        # hidden variable for test
        self.scatter = scatter()

        # spectrometer
        self.R = 300 # scattering radius in mm
        self.MCP_size = np.array([25, 25]) # width x height in mm
        self.MCP_pixel = np.array([128, 128])

    #########################################################################
    ############################## preparation ##############################
    #########################################################################

    def set_peak_1(self, hkl, position):
        [th, tth, z] = position
        self.scatter.set_peak_1(hkl, [tth, 0, self.scatter.beam.beam_energy], [th, 0, 0])
        return

    def set_peak_2(self, hkl, position):
        [th, tth, z] = position
        self.scatter.set_peak_2(hkl, [tth, 0, self.scatter.beam.beam_energy], [th, 0, 0])
        return

    def initialize(self):
        self.IO_initialize()
        self.scatter.cal_UB_matrix()
        self.save_configuration()
        return

    def help(self):
        '''
            Print the help document.
        '''
        print("============ Manual ============")
        print("----- A. Preparations -----")
        print("Assume the object name is <self>, then")
        print("1. Call self.scatter.set_lat([a, b, c], [alpha, beta, gamma]).")
        print("2. Call self.scatter.set_beam(beam_energy).")
        print("3. Call self.set_peak_1(hkl, position) &")
        print("        self.set_peak_2(hkl, position).")
        print("4. Call self.initialize().")
        return

    def info(self):
        '''
            Print current configurations.
        '''
        print("============ Project Name ============")
        print(f"Name: {self.PROJECT_NAME}")
        print()
        self.scatter.info()
        print()
        print("======= Spectrometer Constant ========")
        print(f"scattering radius = {self.R} mm")
        print(f"MCP size = {self.MCP_size} mm x mm")
        print(f"MCP pixels = {self.MCP_pixel}")
        return

    ######################################################################
    ############################## data I/O ##############################
    ######################################################################

    def IO_initialize(self):
        '''
            Initialize input / output system. Create directories.
            IMPORTANT:
                Only call this method after specified self.PROJECT_NAME!
        '''
        self.create_directory()
        return

    def create_directory(self):
        '''
            Check and create the directories.
            Data will be saved in these directories.
        '''
        DIR1 = './' + self.PROJECT_NAME + '/Data/Scans'
        DIR2 = './' + self.PROJECT_NAME + '/Data/Scans_MCP'
        DIR3 = './' + self.PROJECT_NAME + '/Data/Rawdata'
        DIR4 = './' + self.PROJECT_NAME + '/Data/MCP_images'
        DIR5 = './' + self.PROJECT_NAME + '/Data/hklc_data'
        if not os.path.exists(DIR1):
            os.makedirs(DIR1)
        if not os.path.exists(DIR2):
            os.makedirs(DIR2)
        if not os.path.exists(DIR3):
            os.makedirs(DIR3)
        if not os.path.exists(DIR4):
            os.makedirs(DIR4)
        if not os.path.exists(DIR5):
            os.makedirs(DIR5)

        # create log file
        #self.create_log()
        #logging.info('Initialize project.')
        #logging.info(f"Create directories: {DIR1}, {DIR2}, {DIR3}.")
        return

    def split_data(self, INPUT_FILENAME, scan_no_range = 0):
        '''
            Use program "unspec" to split raw data into series of scans.
            scan_no_range = [a, b] : extract scans from #a to #b, optional to avoid the lastest broken scan.
        '''
        INPUT_DIR = './' + self.PROJECT_NAME + '/Data/Rawdata/'
        OUTPUT_DIR = './' + self.PROJECT_NAME + '/Data/Scans/'
        COMMAND = "./unspec " + INPUT_DIR + INPUT_FILENAME + ' ' + OUTPUT_DIR + "scan_ dat -3 -v"

        if scan_no_range != 0:
            COMMAND = COMMAND + " -r " + str(scan_no_range[0]) + ' ' + str(scan_no_range[1])

        os.system(COMMAND)
        print("Scans extracted.")
        return

    def split_data_mcp(self, INPUT_FILENAME):
        '''
            Split mcp raw data into series of scans.
            scan_no_range = [a, b] : extract scans from #a to #b, optional to avoid the lastest broken scan.
        '''
        INPUT_DIR = './' + self.PROJECT_NAME + '/Data/Rawdata/'
        OUTPUT_DIR = './' + self.PROJECT_NAME + '/Data/Scans_MCP/'

        with open(INPUT_DIR + INPUT_FILENAME) as f:
            lines = f.read()

        scans = lines.split("\n\n")

        for scan in scans:
            scan_num = scan[0:10].split(" ")[1] # wicked trick to extract scan No.
            with open(OUTPUT_DIR + 'scan_mcp_' + str(scan_num).zfill(3) + '.dat', 'w') as f:
                if scan[0] == '\n':
                    scan = scan[1:] # remove empty line
                f.write(scan)
                print(f"Extracted scan #{scan_num}.")
        return

    def scan_filename(self, scan_num):
        '''
            Generate filename for scan #scan_num.
        '''
        FILENAME = './' + self.PROJECT_NAME + '/Data/Scans/scan_' + str(scan_num).zfill(3) + '.dat'
        return FILENAME

    def scan_mcp_filename(self, scan_num):
        '''
            Generate filename for scan #scan_num.
        '''
        FILENAME = './' + self.PROJECT_NAME + '/Data/Scans_MCP/scan_mcp_' + str(scan_num).zfill(3) + '.dat'
        return FILENAME

    def scan_mcp_animation_filename(self, scan_num):
        DIR_movie = './' + self.PROJECT_NAME + '/Data/MCP_images/movie_scan_' + str(scan_num).zfill(3) + '.mp4'
        return DIR_movie

    def scan_info(self, scan_num):
        '''
            Print information from scan #scan_num.
        '''
        FILENAME = self.scan_filename(scan_num)

        with open(FILENAME) as f:
            for i in range(5):
                print(f.readline())
        return

    def read_ascan(self, scan_num):
    
        '''
            Read SPEC ascan data from a single data file with specific scan #.
            
            args:
                FILEPATH : [string] to indicate where the data file is located.
                FILENAME : [string] to specific the data file name.
                scan_num : [int] to specific the scan number.
                
            returns:
                [th, I0_BD3, TEY, MCP, pm3] : [list]
                
            example:
                [th, I0_BD3, TEY, MCP, pm3] = spec_scan_reader("./data/", "Sample_A", 56)
        '''
        FILENAME = self.scan_filename(scan_num)

        # read data
        data_origin = np.loadtxt(FILENAME)

        # process data
        data = np.transpose(data_origin)
        th = data[0]
        I0_BD3 = data[6]
        TEY = data[7]
        MCP = data[8]
        pm3 = data[12]
        
        return [th, I0_BD3, TEY, MCP, pm3]

    def read_a2scan(self, scan_num):
    
        '''
            Read SPEC ascan data from a single data file with specific scan #.
            
            args:
                FILEPATH : [string] to indicate where the data file is located.
                FILENAME : [string] to specific the data file name.
                scan_num : [int] to specific the scan number.
                
             returns:
                [tth, th, I0_BD3, TEY, MCP, pm3] : [list]
                
            example:
                [tth, th, I0_BD3, TEY, MCP, pm3] = spec_th2th_reader("./data/", "Sample_A", 56)
        '''
        FILENAME = self.scan_filename(scan_num)

        # read data
        data_origin = np.loadtxt(FILENAME)

        # process data
        data = np.transpose(data_origin)
        tth = data[0]
        th = data[1]
        I0_BD3 = data[7]
        TEY = data[8]
        MCP = data[9]
        pm3 = data[13]
        
        return [tth, th, I0_BD3, TEY, MCP, pm3]

    def save_data_hklc(self, hklc_data, FILENAME):
        '''
            Save hklc_data to file (appending mode).
            Do not attach suffix to FILENAME!
        '''
        DIR = './' + self.PROJECT_NAME + '/Data/hklc_data/'

        with open(DIR + FILENAME + '.csv', 'a', newline='') as csvfile:
            datawriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for data_point in hklc_data:
                datawriter.writerow(data_point)
        return

    def read_data_hklc(self, FILENAME):
        '''
            Read hklc_data from file.
            Do not attach suffix to FILENAME!
        '''
        DIR = './' + self.PROJECT_NAME + '/Data/hklc_data/'

        hklc_data = []
        with open(DIR + FILENAME + '.csv', newline='') as csvfile:
            datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
            line_num = 0
            for row in datareader:
                hklc_data.append([float(row[0]), float(row[1]), float(row[2]), float(row[3])])
                line_num = line_num + 1
        print(f'Processed {line_num} lines.')

        return np.array(hklc_data)

    def save_configuration(self):
        '''
            Save configurations to file.
        '''
        FILENAME = './' + self.PROJECT_NAME +'/configuration.csv'
        
        # write data
        with open(FILENAME, 'w', newline='') as csvfile:
            datawriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            datawriter.writerow(['========== Lattice Constant =========='])
            datawriter.writerow(['lengths', self.scatter.crystal.lattice_lengths])
            datawriter.writerow(['angles', self.scatter.crystal.lattice_angles])
            datawriter.writerow([' '])
            datawriter.writerow(['================ Beam ================'])
            datawriter.writerow(['beam energy', self.scatter.beam.beam_energy])
            datawriter.writerow([' '])
            datawriter.writerow(['============ Bragg Peaks ============='])
            datawriter.writerow(['peak 1'])
            datawriter.writerow(['hkl index', self.scatter.peak_1[0]])
            datawriter.writerow(['theta, phi, energy', self.scatter.peak_1[1]])
            datawriter.writerow(['omega, mu, nu', self.scatter.peak_1[2]])
            datawriter.writerow(['peak 2'])
            datawriter.writerow(['hkl index', self.scatter.peak_2[0]])
            datawriter.writerow(['theta, phi, energy', self.scatter.peak_2[1]])
            datawriter.writerow(['omega, mu, nu', self.scatter.peak_2[2]])
            datawriter.writerow(['UB matrix', self.scatter.UB_matrix])

        return

    #################################################################################
    ############################## lattice calculation ##############################
    #################################################################################

    def angles_2_hkl(self, position_data):
        '''
            Convert position data to hkl index using UB matrix.
        '''
        [th, tth, z] = position_data
        gamma = np.arctan(z / self.R)
        theta = np.rad2deg(np.arccos(np.cos(gamma) * np.cos(np.deg2rad(tth))))
        phi = np.rad2deg(np.arcsin(np.sin(gamma) / np.sin(np.deg2rad(theta))))

        angles = np.array([theta, phi, self.scatter.beam.beam_energy])
        rotations = np.array([th, 0, 0])
        hkl = self.scatter.position2hkl(angles, rotations)

        return hkl

    ################################################################################
    ############################## process MCP data ################################
    ################################################################################

    def img_data_scan_mcp(self, scan_num):
        '''
            Extract counts and position data from MCP images in one scan.
        '''
        FILENAME = self.scan_mcp_filename(scan_num)

        with open(FILENAME) as f:
            datablock = f.read()
        db = datablock.split("#C TwoTheta\t Detz\n")

        command = db[0].split('\n')[0].split(" ")
        if command[2] == 'ascan':
            th_min = float(command[5])
            th_max = float(command[6])
            th_interval = int(command[8])
            th = np.linspace(th_min, th_max, th_interval + 1)

        print("========== Scan Info. ==========")
        print(db[0][0:-1])
        print("================================")

        img_no = len(db[1:]) # number of images
        # row_no = int(img_no / 2 + 0.5)

        imgs_data = np.zeros([img_no, self.MCP_pixel[0], self.MCP_pixel[1]])
        positions_data = np.zeros([img_no, self.MCP_pixel[0], self.MCP_pixel[1], 3])

        img_i = 0
        for snap in db[1:]:
            # processing image data
            img_data_raw = snap.split("\n#@IMG\n")[1].split('\n')
            if len(img_data_raw) == self.MCP_pixel[1] + 1:
                img_data_raw = img_data_raw[0:-1] # ignore the final empty line

            i = 0
            for line in img_data_raw:
                if line[-1] == ' ':
                    line = line[0:-1] # ignore the final space
                imgs_data[img_i][i] = list(map(int, line.split(' ')))
                i = i + 1

            # processing position data
            position_data_raw = snap.split("\n#@IMG\n")[0].split('\n')
            tth = np.zeros([self.MCP_pixel[0]])
            z = np.zeros([self.MCP_pixel[1]])
            i = 0
            for position in position_data_raw:
                [tth[i], z[i]] = list(map(float, position.split(' ')))
                i = i + 1

            for i in range(self.MCP_pixel[0]):
                for j in range(self.MCP_pixel[1]):
                    positions_data[img_i][i][j] = [th[img_i], tth[i], z[j]]

            img_i = img_i + 1

        return (imgs_data, positions_data)

    def hkl_data_scan_mcp(self, scan_num):
        '''
            Extract counts and position data from MCP images in one scan,
            and then map the postition data to hkl space.
        '''
        print("Reading MCP data ...")
        (imgs_data, positions_data) = self.img_data_scan_mcp(scan_num)
        print("Done.")
        print("Converting MCP data to HKL space ...")
        img_no = len(imgs_data)
        hkl_positions_data = np.zeros([img_no, self.MCP_pixel[0], self.MCP_pixel[1], 3])
        for img_i in range(img_no):
            for i in range(self.MCP_pixel[0]):
                for j in range(self.MCP_pixel[1]):
                    hkl_positions_data[img_i][i][j] = self.angles_2_hkl(positions_data[img_i][i][j])
        print("Done.")
        return (imgs_data, hkl_positions_data)

    def hklc_data_combiner_mcp(self, imgs_data, hkl_positions_data):
        '''
            Prepare data for scatter in 3D hkl space.
        '''
        # calculate data amounts
        N = len(imgs_data) * self.MCP_pixel[0] * self.MCP_pixel[1]

        # create empty arrays. hklc_data = [h, k, l, counts]
        hklc_data = np.zeros([N, 4])

        k = 0
        for img_i in range(len(imgs_data)):
            for i in range(self.MCP_pixel[0]):
                for j in range(self.MCP_pixel[1]):
                    hklc_data[k] = np.append(hkl_positions_data[img_i][i][j], imgs_data[img_i][i][j])
                    k = k + 1
        return hklc_data

    def hklc_data_filter(self, hklc_data, threshold):
        '''
            Ignore data with low counts.
        '''
        hklc_data_filtered = []

        i = 0
        for i in range(len(hklc_data)):
            if hklc_data[i][3] >= threshold:
                hklc_data_filtered.append(hklc_data[i])
        return np.array(hklc_data_filtered)

    def hklc_data_scan_mcp(self, scan_num):
        '''
            Prepare [h, k, l, counts] data for 3d scatter plot.
        '''
        start_time = time.time()
        # get data from single scan file
        (imgs_data, hkl_positions_data) = self.hkl_data_scan_mcp(scan_num)

        # processing and combine data
        hklc_data = self.hklc_data_combiner_mcp(imgs_data, hkl_positions_data)

        elapsed_time = time.time() - start_time
        print("Time consuming: {0:.3f}s.".format(elapsed_time))

        return hklc_data

    #############################################################################
    ############################## MCP visualize ################################
    #############################################################################

    def display_scan_mcp(self, scan_num, v_range, column_no, show, save):
        '''
            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
            v_max = -2: use different colorbar for each image
        '''
        # read data
        (imgs_data, positions_data) = self.img_data_scan_mcp(scan_num)
        [v_min, v_max] = v_range

        img_no = len(imgs_data) # number of images
        row_no = int(img_no / column_no) + (img_no % column_no > 0) # round up
        
        fig = plt.figure(figsize=(20, 20 / column_no * row_no))

        if v_max == -1: # use highest count as vmax
            v_max = np.amax(imgs_data)
        elif v_max == -2: # use different colorbar for each image
            v_max = None

        for i in range(img_no):
            ax = fig.add_subplot(row_no, column_no, i + 1)
            ax.set_title('snap ' + str(i))
            img = plt.imshow(imgs_data[i].T, origin='lower', vmin=v_min, vmax=v_max)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            plt.colorbar(img, cax)

        fig.tight_layout()

        plt.show()

        return

    def display_scan_mcp_1(self, scan_num, v_range, snap_no, fig_size, show, save):
        '''
            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
            v_max = -2: use different colorbar for each image
        '''
        # read data
        (imgs_data, positions_data) = self.img_data_scan_mcp(scan_num)
        [v_min, v_max] = v_range
        
        fig = plt.figure(figsize=fig_size)

        if v_max == -1: # use highest count as vmax
            v_max = np.amax(imgs_data)
        elif v_max == -2: # use different colorbar for each image
            v_max = None

        ax = fig.add_subplot(111)
        ax.set_title('scan ' + str(scan_num) + ' snap ' + str(snap_no))
        img = plt.imshow(imgs_data[snap_no].T, origin='lower', vmin=v_min, vmax=v_max)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(img, cax)

        fig.tight_layout()

        if save == 1:
            DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/scan_' + str(scan_num).zfill(3) + '_' + str(snap_no).zfill(2) + '.png'
            plt.savefig(DIR, dpi = 150, format = 'png')

        if show == 1:
            plt.show()
        elif show == 0:
            plt.close(fig)

        return

    def animation_scan_mcp(self, scan_num, v_range = [0, -1], fig_size = (10, 10), font_size = 20, show = 1, clean = 1):
        '''
            Generate animations with specified MCP scan_num.
            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
        '''
        start_time = time.time()

        # read data
        (imgs_data, positions_data) = self.img_data_scan_mcp(scan_num)
        [v_min, v_max] = v_range
        if v_max == -1: # use highest count as vmax
            v_max = np.amax(imgs_data)

        # prepare directory
        DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/scan_' + str(scan_num).zfill(3)
        if not os.path.exists(DIR):
            os.makedirs(DIR)
        DIR_movie = './' + self.PROJECT_NAME + '/Data/MCP_images/movie_scan_' + str(scan_num).zfill(3) + '.mp4'

        # plotting
        print("Start generating images...")
        for snap_no in range(len(imgs_data)):
            fig = plt.figure(figsize=fig_size)
            plt.rcParams.update({'font.size': font_size})

            ax = fig.add_subplot(111)
            ax.set_title('snap ' + str(snap_no))
            img = plt.imshow(imgs_data[snap_no], vmin=v_min, vmax=v_max)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            plt.colorbar(img, cax)

            fig.tight_layout()
            plt.savefig(DIR + '/snap_' + str(snap_no).zfill(2) + '.png', dpi = 150, format = 'png')
            plt.close()
        print("Done.")

        # generate animation using ffmpeg
        print("Start generating animation...")
        os.system("ffmpeg -r 5 -i " + DIR + "/snap_%02d.png" 
                  + " -vcodec mpeg4 -y " + DIR_movie)
        print("Done.")
        # delete all figures
        if clean == 1:
            shutil.rmtree(DIR)
            print("All images are deleted.")

        elapsed_time = time.time() - start_time
        print("Time consuming: {0:.3f}s.".format(elapsed_time))

        if show == 1:
            return Video(DIR_movie, width=600, height=600)

        return

    def new_colormaps(self):
        '''
            Define several new colormaps
        '''
        from matplotlib.colors import LinearSegmentedColormap

        ########## jet_alpha ##########
        # get colormap
        ncolors = 256
        color_array = plt.get_cmap('jet')(range(ncolors))

        # change alpha values
        color_array[:,-1] = np.linspace(0,1,ncolors)

        # create a colormap object
        map_object = LinearSegmentedColormap.from_list(name='jet_alpha',colors=color_array)

        # register this new colormap with matplotlib
        plt.register_cmap(cmap=map_object)

        ########## jet_alpha2 ##########
        # get colormap
        ncolors = 256
        color_array = plt.get_cmap('jet')(range(ncolors))

        # change alpha values
        a = np.linspace(0,1,ncolors)
        alpha_array = 0.5 * (1 - np.cos(np.pi * a))
        color_array[:,-1] = alpha_array

        # create a colormap object
        map_object = LinearSegmentedColormap.from_list(name='jet_alpha2',colors=color_array)

        # register this new colormap with matplotlib
        plt.register_cmap(cmap=map_object)

        return

    def scatter_HKL(self, hklc_data, threshold, log=0, marker_size=15, colormap='jet_alpha', equal_length=0):
        '''
            Generating 3D scatter plot in HKL space.
        '''
        # register new colormap
        self.new_colormaps()

        # filter data
        hklc_data_filtered = self.hklc_data_filter(hklc_data, threshold)

        # prepare
        db = np.transpose(hklc_data_filtered)

        fig = plt.figure(figsize = (10, 7))
        ax = fig.add_subplot(111, projection='3d')
        if log == 0:
            img = ax.scatter(db[0], db[1], db[2], c=db[3], cmap = colormap)
        elif log == 1:
            img = ax.scatter(db[0], db[1], db[2], c=np.log(db[3]), cmap = colormap, edgecolor = 'none', s=marker_size)
        fig.colorbar(img, shrink=0.85)
        ax.set_xlabel('H')
        ax.set_ylabel('K')
        ax.set_zlabel('L')

        if equal_length == 1:
            [X, Y, Z] = db[0:3]
            max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0

            mid_x = (X.max()+X.min()) * 0.5
            mid_y = (Y.max()+Y.min()) * 0.5
            mid_z = (Z.max()+Z.min()) * 0.5
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)

        fig.tight_layout()
        #plt.savefig("HKL_space.png", dpi = 150, format = 'png')
        plt.show()
        return

    #######################################################################
    ############################## visualize ##############################
    #######################################################################

    def plot_ascan(self, scan_nums, DETECTOR, I0_BD3, fig_size, font_size, data_legends, LEGEND_PREFIX, LEGEND_SUFFIX, TITLE = ''):
        '''
            Plot a set of processed SPEC data together with legends. 
            
            Author: Wenjie Chen 
            E-mail: wenjiechen@pku.edu.cn
            
            args:
                datablock : [list] data read from SPEC data file with function "spec_data_reader".
                SCAN_FORMAT : [string] indicate the scan format, "ascan" or "a2scan".
                VARIABLE : [string] indicate the VARIABLE, "th" or "tth".
                DETECTOR : [string] must be chosen from "TEY", "MCP", "pm3".
                I0_BD3 : [int] 0 - do not divided by I0_BD3, 1 - divided by I0_BD3.
                data_legends : [list] a set of legends, [20, 40, 60], e.g.
                LEGEND_PREFIX : [string] prefix to data_legends.
                LEGEND_SUFFIX : [string] suffix to data_legends.
                TITLE : [string] figure title.
                
            returns:
                a figure with multiple curves marked with legends
                
            example:
                spec_plot(datablock, "a2scan", "tth", "MCP", 1, [100, 101, 102], "scan #", "", "differnet scans")
                spec_plot(datablock, "ascan", "th", "MCP", 0, [20, 50, 80], "T = ", " K", "Temperature dependence")
        ''' 
        
        plt.figure(figsize=fig_size)
        plt.rcParams.update({'font.size': font_size})

        datablock = []
        for scan_num in scan_nums:
            datablock.append(self.read_ascan(scan_num))
            
        if DETECTOR == "TEY":
            det = 2
        elif DETECTOR == "MCP":
            det = 3
        elif DETECTOR == "pm3":
            det = 4
        else:
            raise ValueError('DETECTOR must be "TEY", "MCP" or "pm3"!')
                
        if I0_BD3 == 0:
            i = 0
            for legend in data_legends:
                plt.plot(datablock[i][0], datablock[i][det], label = LEGEND_PREFIX + str(legend) + LEGEND_SUFFIX)
                i = i + 1
        elif I0_BD3 == 1:
            i = 0
            for legend in data_legends:
                plt.plot(datablock[i][0], datablock[i][det]/datablock[i][1], label = LEGEND_PREFIX + str(legend) + LEGEND_SUFFIX)
                i = i + 1
        else:
            raise ValueError('IO_BD3 must be 0 or 1!')

        plt.xlabel('th')
        plt.ylabel(DETECTOR)
        plt.legend()
        plt.title(TITLE)
        plt.show()
        return

    def plot_a2scan(self, scan_nums, DETECTOR, I0_BD3, fig_size, font_size, data_legends, LEGEND_PREFIX, LEGEND_SUFFIX, TITLE = ''):
        '''
            Plot a set of processed SPEC data together with legends. 
            
            Author: Wenjie Chen 
            E-mail: wenjiechen@pku.edu.cn
            
            args:
                datablock : [list] data read from SPEC data file with function "spec_data_reader".
                SCAN_FORMAT : [string] indicate the scan format, "ascan" or "a2scan".
                VARIABLE : [string] indicate the VARIABLE, "th" or "tth".
                DETECTOR : [string] must be chosen from "TEY", "MCP", "pm3".
                I0_BD3 : [int] 0 - do not divided by I0_BD3, 1 - divided by I0_BD3.
                data_legends : [list] a set of legends, [20, 40, 60], e.g.
                LEGEND_PREFIX : [string] prefix to data_legends.
                LEGEND_SUFFIX : [string] suffix to data_legends.
                TITLE : [string] figure title.
                
            returns:
                a figure with multiple curves marked with legends
                
            example:
                spec_plot(datablock, "a2scan", "tth", "MCP", 1, [100, 101, 102], "scan #", "", "differnet scans")
                spec_plot(datablock, "ascan", "th", "MCP", 0, [20, 50, 80], "T = ", " K", "Temperature dependence")
        ''' 
        
        plt.figure(figsize=fig_size)
        plt.rcParams.update({'font.size': font_size})

        datablock = []
        for scan_num in scan_nums:
            datablock.append(self.read_a2scan(scan_num))
            
        if DETECTOR == "TEY":
            det = 2
        elif DETECTOR == "MCP":
            det = 3
        elif DETECTOR == "pm3":
            det = 4
        else:
            raise ValueError('DETECTOR must be "TEY", "MCP" or "pm3"!')
                
        if I0_BD3 == 0:
            i = 0
            for legend in data_legends:
                plt.plot(datablock[i][0], datablock[i][det], label = LEGEND_PREFIX + str(legend) + LEGEND_SUFFIX)
                i = i + 1
        elif I0_BD3 == 1:
            i = 0
            for legend in data_legends:
                plt.plot(datablock[i][0], datablock[i][det]/datablock[i][2], label = LEGEND_PREFIX + str(legend) + LEGEND_SUFFIX)
                i = i + 1
        else:
            raise ValueError('IO_BD3 must be 0 or 1!')

        plt.xlabel('tth')
        plt.ylabel(DETECTOR)
        plt.legend()
        plt.title(TITLE)
        plt.show()
        return
