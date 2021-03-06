########## INFO ##########
print("########################################")
print("Project: spec_data_analyser")
print("Version: 1.4.5 - +[SDD_ROI sum 2 average]")
print("Last Update: 2021.05.25")
print("----------------------------------------")
print("Author: Wenjie Chen")
print("E-mail: wenjiechen@pku.edu.cn")
print("########################################")
##########################

import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import random

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

        # MCP detector
        self.MCP_size = np.array([25, 25]) # width x height in mm
        self.MCP_pixel = np.array([128, 128])
        self.MCP_response_img_data = np.zeros([128, 128])
        self.MCP_effective_pixels = 128 * 128
        self.MCP_nonlinear_repair_method_dict = {'NONE' : self.img_data_mcp_nonlinear_NONE,
                                                 '2nd' : self.img_data_mcp_nonlinear_2nd,
                                                 '3rd' : self.img_data_mcp_nonlinear_3rd}
        self.MCP_nonlinear_repair_method_key = 'NONE'
        self.MCP_nonlinear_repair_method_params = None

        # global switch
        self.ctps_key = False # all relative variables to 'per second' unit
        self.I0_BD3_key = False # all relative variables divided by I0_BD3
        self.MCP_nonlinear_repair_key = False
        self.MCP_non_uniform_repair_key = False

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
        print()
        print("======= MCP Detector Constant ========")
        print(f"MCP size = {self.MCP_size} mm x mm")
        print(f"MCP pixels = {self.MCP_pixel}")
        print()
        print("========= Global Switch List =========")
        print(f"Use per second unit: {self.ctps_key}")
        print(f"Divided by I0_BD3: {self.I0_BD3_key}")
        print(f"Repair MCP nonlinear response: {self.MCP_nonlinear_repair_key}")
        print(f"Repair MCP non-uniform response: {self.MCP_non_uniform_repair_key}")
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
        DIR6 = './' + self.PROJECT_NAME + '/Data/MCP_data'
        DIR7 = './' + self.PROJECT_NAME + '/Data/Scans_SDD'

        DIRs = [DIR1, DIR2, DIR3, DIR4, DIR5, DIR6, DIR7]

        for DIR in DIRs:
            if not os.path.exists(DIR):
                os.makedirs(DIR)
        # if not os.path.exists(DIR1):
        #     os.makedirs(DIR1)
        # if not os.path.exists(DIR2):
        #     os.makedirs(DIR2)
        # if not os.path.exists(DIR3):
        #     os.makedirs(DIR3)
        # if not os.path.exists(DIR4):
        #     os.makedirs(DIR4)
        # if not os.path.exists(DIR5):
        #     os.makedirs(DIR5)
        # if not os.path.exists(DIR6):
        #     os.makedirs(DIR6)

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

    def split_data_mcp(self, INPUT_FILENAME, hide_details = 1):
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
                if hide_details != 1:
                    print(f"Extracted scan #{scan_num}.")
        print("Scans MCP data extracted.")
        return

    def split_data_sdd(self, INPUT_FILENAME, hide_details = 1):
        '''
            Split sdd raw data into series of scans.
            scan_no_range = [a, b] : extract scans from #a to #b, optional to avoid the lastest broken scan.
        '''
        INPUT_DIR = './' + self.PROJECT_NAME + '/Data/Rawdata/'
        OUTPUT_DIR = './' + self.PROJECT_NAME + '/Data/Scans_SDD/'
        
        with open(INPUT_DIR + INPUT_FILENAME) as f:
            lines = f.read()
            
        scans = lines.split("\n\n")
        
        for scan in scans:
            scan_num = scan[0:10].split(" ")[1] # wicked trick to extract scan No.
            with open(OUTPUT_DIR + 'scan_sdd_' + str(scan_num).zfill(3) + '.dat', 'w') as f:
                if scan[0] == '\n':
                    scan = scan[1:] # remove empty line
                f.write(scan)
                if hide_details != 1:
                    print(f"Extracted scan #{scan_num}.")
        print("Scans SDD data extracted.")

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

    def scan_sdd_filename(self, scan_num):
        '''
            Generate filename for scan #scan_num.
        '''
        FILENAME = './' + self.PROJECT_NAME + '/Data/Scans_SDD/scan_sdd_' + str(scan_num).zfill(3) + '.dat'
        return FILENAME

    def read_scan(self, scan_num):
        '''
            Read data, return column name and data block.
        '''
        FILENAME = self.scan_filename(scan_num)

        # read data
        with open(FILENAME) as f:
            raw_data = f.readlines()
            scan_info = raw_data[:5]
            column_name = raw_data[5].split()[1:]
            data_block_raw = raw_data[6:]

        # process data
        data_num = []
        for line in data_block_raw:
            data_txt = line.split(' ')
            data_txt[-1] = data_txt[-1][:-1]
            data_num.append(list(map(float, data_txt)))
            
        data_block = np.transpose(np.array(data_num))
        return [column_name, data_block]

    def scan_value(self, scan_num, variable):
        '''
            Return the data array named <variable> in scan #scan_num.
        '''
        [column_name, data_block] = self.read_scan(scan_num)
        return data_block[column_name.index(variable)]

    def read_ascan(self, scan_num):
    
        '''
            Read SPEC ascan data from a single data file with specific scan #.
            Shall be deleted in the future.
            
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
        MCP = data[-1]
        pm3 = data[12]

        H = data[1]
        K = data[2]
        L = data[3]
        
        return [th, I0_BD3, TEY, MCP, pm3, H, K, L]
        
    def read_a2scan(self, scan_num):
    
        '''
            Read SPEC ascan data from a single data file with specific scan #.
            Shall be deleted in the future.
            
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

        H = data[1]
        K = data[2]
        L = data[3]
        
        return [tth, th, I0_BD3, TEY, MCP, pm3, H, K, L]

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

    ###############################################################################
    ############################## extract scan info ##############################
    ###############################################################################

    def T_info(self, scan_nums):
        '''
            Display temperature information for scans.
        '''
        for scan_num in scan_nums:
            T = np.mean(self.scan_value(scan_num, 'Sample_T'))
            print(f"Scan #{scan_num}, T = {T}K")
        return

    def column_name_info(self, scan_num):
        '''
            Display column name information for scan.
        '''
        [column_name, data_block] = self.read_scan(scan_num)
        print(column_name)
        return

    def command_info(self, scan_nums):
        '''
            Display scan command for scans.
        '''
        for scan_num in scan_nums:
            FILENAME = self.scan_filename(scan_num)

            with open(FILENAME) as f:
                COMMAND = f.readline()[0:-1]
                T = np.mean(self.scan_value(scan_num, 'Sample_T'))
                T = round(T, 1)
                print(f"{COMMAND}, T_mean = {T}K")
        return

    def scan_info(self, scan_num):
        '''
            Print information from scan #scan_num.
        '''
        FILENAME = self.scan_filename(scan_num)

        with open(FILENAME) as f:
            for i in range(5):
                print(f.readline())
        return

    def get_exposure_time(self, scan_num):
        FILENAME = self.scan_filename(scan_num)
        with open(FILENAME) as f:
            COMMAND = f.readline()[0:-1]
            strings = COMMAND.split(' ')
            if strings[-1] == 'fixQ': # special care for EfixQ scan
                exposure_time = float(strings[-2])
            elif strings[4] == 'tscan': # special care for tscan
                exposure_time = float(strings[-2])
            else: # for most scans
                exposure_time = float(strings[-1])
        return exposure_time

    def get_T(self, scan_nums):
        T = []
        for scan_num in scan_nums:
            T.append(np.mean(self.scan_value(scan_num, 'Sample_T')))
        return np.array(T)

    def get_pol_angle(self, scan_nums):
        pol = []
        for scan_num in scan_nums:
            FILENAME = self.scan_filename(scan_num)
            with open(FILENAME) as f:
                for i in range(5):
                    line = f.readline()
            pol.append(float(line.split(' ')[4]))
        return np.array(pol)

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

    #################################################################################
    ############################## process scan data ################################
    #################################################################################

    def average_scan_value(self, scan_nums, variable, I0 = 1):
        flag = 0
        if I0 == 1: # divided by I0_BD3
            for scan_num in scan_nums:
                data = self.scan_value(scan_num, variable)
                I0_BD3 = self.scan_value(scan_num, 'I0_BD3')
                if flag == 0:
                    data_average = np.zeros([len(data)])
                    flag = 1
                data_average = data_average + data / I0_BD3
            data_average = data_average / len(scan_nums)
        else: # not divided by I0_BD3
            for scan_num in scan_nums:
                data = self.scan_value(scan_num, variable)
                if flag == 0:
                    data_average = np.zeros([len(data)])
                    flag = 1
                data_average = data_average + data
            data_average = data_average / len(scan_nums)
        return data_average

    def average_ascan(self, scan_nums):
        '''
            Shall be deleted in the future.
        '''
        flag = 0
        for scan_num in scan_nums:
            [th, I0_BD3, TEY, MCP, pm3, H, K, L] = self.read_ascan(scan_num)
            if flag == 0:
                TEY_average = np.zeros([len(TEY)])
                MCP_average = np.zeros([len(MCP)])
                pm3_average = np.zeros([len(pm3)])
                flag = 1
            TEY_average = TEY_average + TEY / I0_BD3
            MCP_average = MCP_average + MCP / I0_BD3
            pm3_average = pm3_average + pm3 / I0_BD3
        TEY_average = TEY_average / len(scan_nums)
        MCP_average = MCP_average / len(scan_nums)
        pm3_average = pm3_average / len(scan_nums)
        return [th, TEY_average, MCP_average, pm3_average, H, K, L]


    ################################################################################
    ############################## process MCP data ################################
    ################################################################################

    def img_data_mcp_nonlinear_NONE(self, img_data, params):
        return img_data

    def img_data_mcp_nonlinear_2nd(self, img_data, params):
        a = params
        MCP_total = np.sum(img_data)
        MCP_total_repaired = -0.5 / a * (1 - np.sqrt(1 + 4 * a * MCP_total))
        scale = MCP_total_repaired / MCP_total
        img_data_repaired = img_data * scale
        return img_data_repaired

    def img_data_mcp_nonlinear_3rd(self, img_data, params):
        print('Debug: 3rd order repairing method has been called.')
        img_data_repaired = img_data
        return img_data_repaired

    def img_data_mcp_nonlinear_repair(self, img_data, METHOD, params):
        '''
            Repair MCP nonlinear response with METHOD from
                '2nd', '3rd' and '???'

            !!! Input img_data should have the unit of counts per second! !!!

            Return: repaired img_data
        '''
        img_data_repaired = self.MCP_nonlinear_repair_method_dict[METHOD](img_data, params)
        return img_data_repaired

    def img_data_round_filter(self, img_data, center, radius):
        '''
            Apply an artificial hole-filter (to get rid of noise outside of the detecting area).
        '''
        img_data_filtered = img_data.copy()
        R_square = radius ** 2 # accelerate speed
        X, Y = center # for CLS, center = [63, 66], radius = 58
        for i in range(self.MCP_pixel[0]):
            for j in range(self.MCP_pixel[0]):
                if (i - X) ** 2 + (j - Y) ** 2 > R_square:
                    img_data_filtered[i][j] = 0
        return img_data_filtered

    def img_data_save(self, img_data, FILENAME):
        DIR = './' + self.PROJECT_NAME + '/Data/MCP_data/' + FILENAME
        np.save(DIR, img_data)
        print("MCP image data has been successfully saved.")
        return

    def img_data_load_mcp_response(self, FILENAME, rescale = True):
        '''
            Load MCP response data for correction.
        '''
        DIR = './' + self.PROJECT_NAME + '/Data/MCP_data/' + FILENAME
        try:
            self.MCP_response_img_data = np.load(DIR)
            if rescale:
                self.MCP_response_img_data = self.MCP_response_img_data / np.max(self.MCP_response_img_data)
            self.MCP_effective_pixels = len(self.MCP_response_img_data[self.MCP_response_img_data != 0])
            print('MCP response data has been successfully loaded.')
            print(f'MCP effective pixels number is {self.MCP_effective_pixels}.')
        except FileNotFoundError:
            print('File does not exist!')
        return

    def img_scan_mcp(self, scan_num, repair = True):
        '''
            Just extract img data, does not care about position info.
        '''
        FILENAME = self.scan_mcp_filename(scan_num)

        with open(FILENAME) as f:
            datablock = f.read()
        db = datablock.split("#C TwoTheta\t Detz\n")
        
        #print("========== Scan Info. ==========")
        #print(db[0][0:-1])
        #print("================================")

        img_no = len(db[1:]) # number of images
        # row_no = int(img_no / 2 + 0.5)

        imgs_data = np.zeros([img_no, self.MCP_pixel[0], self.MCP_pixel[1]])

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

            # adjust flip
            imgs_data[img_i] = np.flip(imgs_data[img_i], 1)

            img_i = img_i + 1

        # transform the unit to counts per second
        if self.ctps_key:
            imgs_data = imgs_data / self.get_exposure_time(scan_num)
            #imgs_data = imgs_data / 2.0 # debug

        # repair MCP nonlinear response
        if self.MCP_nonlinear_repair_key:
            for i in range(len(imgs_data)):
                imgs_data[i] = self.img_data_mcp_nonlinear_repair(imgs_data[i], self.MCP_nonlinear_repair_method_key, self.MCP_nonlinear_repair_method_params)

        # repair MCP by response data (self.MCP_response_img_data)
        if self.MCP_non_uniform_repair_key and repair and self.MCP_response_img_data.any() != 0: 
            imgs_data = np.true_divide(imgs_data.copy(), self.MCP_response_img_data, where=(self.MCP_response_img_data!=0))
            # I still don't understand why .copy() is needed here. 
            # Without this copy, when self.stps_key = True, "where=(self.MCP_response_img_data!=0)" will not perform correctly.

        # divided by I0_BD3
        if self.I0_BD3_key:
            I0_BD3 = self.scan_value(scan_num, 'I0_BD3') / self.get_exposure_time(scan_num)
            imgs_data = imgs_data / I0_BD3[:, np.newaxis, np.newaxis]

        return imgs_data

    def img_data_scan_mcp(self, scan_num, repair = True):
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
        elif command[2] == 'a2scan':
            th_min = float(command[9])
            th_max = float(command[10])
            th_interval = int(command[12])
            th = np.linspace(th_min, th_max, th_interval + 1)
        elif command[2] == 'hklscan':
            th = self.scan_value(scan_num, 'Theta')
        else:
            th = self.scan_value(scan_num, 'H')
            # th = self.scan_value(scan_num, 'Theta') # should be appliable to all commands

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

            # adjust flip
            imgs_data[img_i] = np.flip(imgs_data[img_i], 1)

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

        # transform the unit to counts per second
        if self.ctps_key:
            imgs_data = imgs_data / self.get_exposure_time(scan_num)
            #imgs_data = imgs_data / 2.0 # debug

        # repair MCP nonlinear response
        if self.MCP_nonlinear_repair_key:
            for i in range(len(imgs_data)):
                imgs_data[i] = self.img_data_mcp_nonlinear_repair(imgs_data[i], self.MCP_nonlinear_repair_method_key, self.MCP_nonlinear_repair_method_params)

        # repair MCP by response data (self.MCP_response_img_data)
        if self.MCP_non_uniform_repair_key and repair and self.MCP_response_img_data.any() != 0: 
            imgs_data = np.true_divide(imgs_data.copy(), self.MCP_response_img_data, where=(self.MCP_response_img_data!=0))
            # I still don't understand why .copy() is needed here. 
            # Without this copy, when self.stps_key = True, "where=(self.MCP_response_img_data!=0)" will not perform correctly.

        # divided by I0_BD3
        if self.I0_BD3_key:
            I0_BD3 = self.scan_value(scan_num, 'I0_BD3') / self.get_exposure_time(scan_num)
            imgs_data = imgs_data / I0_BD3[:, np.newaxis, np.newaxis]

        return (imgs_data, positions_data)

    def img_data_subtract(self, imgs_data_1, imgs_data_2):
        imgs_data_sub = []
        img_no = len(imgs_data_1)
        for img_i in range(img_no):
            imgs_data_sub.append(imgs_data_1[img_i] - imgs_data_2[img_i])
        return np.array(imgs_data_sub)

    def img_data_subtract_bg(self, imgs_data):
        '''
            Subtract counts from the first image.
        '''
        imgs_data_sub_bg = []
        img_no = len(imgs_data)
        for img_i in range(img_no):
            imgs_data_sub_bg.append(imgs_data[img_i] - imgs_data[0])
        return np.array(imgs_data_sub_bg)

    def img_data_ROI_static(self, imgs_data, ROI_range):
        img_no = len(imgs_data)
        MCP_ROI = []
        for img_i in range(img_no):
            count = 0
            for i in range(ROI_range[0][0], ROI_range[0][1]):
                for j in range(ROI_range[1][0], ROI_range[1][1]):
                    count = count + imgs_data[img_i][i][j]
            MCP_ROI.append(count)
        return np.array(MCP_ROI)

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

    def hkl_data_mcp_imgs(self, imgs_data, positions_data):
        '''
            Extract counts and position data from MCP images in one scan,
            and then map the postition data to hkl space.
        '''
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

        # extract counts array
        c = np.transpose(hklc_data)[3]

        # filtering
        hklc_data_filtered = hklc_data[np.where(c >= threshold)]

        return hklc_data_filtered

    def hklc_data_filter_hkl(self, hklc_data, h_range, k_range, l_range):
        '''
            Focus on data in ROI.
        '''
        hklc_data_filtered = []

        i = 0
        for i in range(len(hklc_data)):
            [H, K, L] = hklc_data[i][0:3]
            if H >= h_range[0] and H <= h_range[1]:
                if K >= k_range[0] and K <= k_range[1]:
                    if L >= l_range[0] and L <= l_range[1]:
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

    def hklc_data_mcp_imgs(self, imgs_data, positions_data):
        '''
            Prepare [h, k, l, counts] data for 3d scatter plot.
        '''
        start_time = time.time()

        (imgs_data, hkl_positions_data) = self.hkl_data_mcp_imgs(imgs_data, positions_data)

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
        #(imgs_data, positions_data) = self.img_data_scan_mcp(scan_num)
        imgs_data = self.img_scan_mcp(scan_num)
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

    def display_scan_mcp_1(self, scan_num, snap_no, VARIABLE = '', v_range = [0, -2], sub_bg = 0, fig_size = (8, 8), show = 1, save = 0):
        '''
            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
            v_max = -2: use different colorbar for each image
        '''
        # read data
        #(imgs_data, positions_data) = self.img_data_scan_mcp(scan_num)
        imgs_data = self.img_scan_mcp(scan_num)
        [v_min, v_max] = v_range
        TITLE = 'scan ' + str(scan_num) + ' snap ' + str(snap_no)

        if VARIABLE != '':
            var_data = self.scan_value(scan_num, VARIABLE)
            TITLE = TITLE + ', ' + VARIABLE + ' = ' + str(round(var_data[snap_no], 2))
        
        fig = plt.figure(figsize=fig_size)

        if v_max == -1: # use highest count as vmax
            v_max = np.amax(imgs_data)
        elif v_max == -2: # use different colorbar for each image
            v_max = None

        if sub_bg == 1:
            img_data = imgs_data[snap_no] - imgs_data[0]
        else:
            img_data = imgs_data[snap_no]

        ax = fig.add_subplot(111)
        ax.set_title(TITLE)
        img = plt.imshow(img_data.T, origin='lower', vmin=v_min, vmax=v_max)
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

    def display_mcp(self, img_data, fig_size = (5, 5), font_size = 8, v_range = [0, -1], save = False, FILENAME = 'untitled'):
        '''
            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
            v_max = -2: use different colorbar for each image
        '''
        [v_min, v_max] = v_range
        if v_max == -1: # use highest count as vmax
            v_max = np.amax(img_data)

        fig = plt.figure(figsize=fig_size, dpi = 200)
        plt.rcParams.update({'font.size': font_size})
        ax = fig.add_subplot(111)
        img = plt.imshow(img_data.T, origin='lower', vmin=v_min, vmax=v_max)
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        clb = plt.colorbar(img, cax)
        plt.colorbar(img, cax)
        clb.set_label('intensity (arb. unit)', labelpad=15, rotation=-90)
        fig.tight_layout()

        if save:
            DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/MCP_img_' + FILENAME + '.png'
            plt.savefig(DIR, dpi = 150, format = 'png')

        plt.show()
        return

    def mark_ROI(self, img_data, ROIs, ROI_colors = [], TITLE = '', line_width = 2, font_size = 12, v_range = [0, -1], fig_size = (8, 8), save = 0):
        '''
            Mark ROIs on one MCP image.

            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
            v_max = -2: use different colorbar for each image
        '''
        [v_min, v_max] = v_range
        if v_max == -1: # use highest count as vmax
            v_max = np.amax(img_data)

        fig = plt.figure(figsize=fig_size, dpi = 200)
        ax = fig.add_subplot(111)
        ax.set_title(TITLE)
        img = plt.imshow(img_data.T, origin='lower', vmin=v_min, vmax=v_max)

        if ROI_colors == []:
            ROI_colors = ['r' for i in range(len(ROIs))]

        for i in range(len(ROIs)):
            ROI_range = ROIs[i]
            color = ROI_colors[i]
            line_1 = [[ROI_range[0][0], ROI_range[0][1]], [ROI_range[1][0], ROI_range[1][0]]]
            line_2 = [[ROI_range[0][1], ROI_range[0][1]], [ROI_range[1][0], ROI_range[1][1]]]
            line_3 = [[ROI_range[0][0], ROI_range[0][1]], [ROI_range[1][1], ROI_range[1][1]]]
            line_4 = [[ROI_range[0][0], ROI_range[0][0]], [ROI_range[1][0], ROI_range[1][1]]]
            plt.plot(line_1[0], line_1[1], c=color, linewidth=line_width)
            plt.plot(line_2[0], line_2[1], c=color, linewidth=line_width)
            plt.plot(line_3[0], line_3[1], c=color, linewidth=line_width)
            plt.plot(line_4[0], line_4[1], c=color, linewidth=line_width)
            plt.text(ROI_range[0][0], ROI_range[1][1] + 1, 'ROI ' + str(i+1), color=color, size=font_size)

        plt.xlabel('pixel')
        plt.ylabel('pixel')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        #plt.colorbar(img, cax)
        clb = plt.colorbar(img, cax)
        clb.set_label('intensity (arb. unit)', labelpad=15, rotation=-90)
        fig.tight_layout()
        if save == 1:
            if TITLE == '':
                TITLE = 'untitled_MCP_img'
            DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/' + TITLE + '_ROI' + '.png'
            plt.savefig(DIR, dpi = 300, format = 'png')
        plt.show()
        return

    def animation_scan_mcp(self, scan_num, VARIABLE = '', v_range = [0, -1], sub_bg = 0, fig_size = (10, 10), font_size = 20, show = 1, clean = 1):
        '''
            Generate animations with specified MCP scan_num.
            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
        '''
        start_time = time.time()

        # read data
        #(imgs_data, positions_data) = self.img_data_scan_mcp(scan_num)
        imgs_data = self.img_scan_mcp(scan_num)
        [v_min, v_max] = v_range
        if v_max == -1: # use highest count as vmax
            v_max = np.amax(imgs_data)

        if VARIABLE != '':
            var_data = self.scan_value(scan_num, VARIABLE)

        # prepare directory
        DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/scan_' + str(scan_num).zfill(3)
        if not os.path.exists(DIR):
            os.makedirs(DIR)
        DIR_movie = './' + self.PROJECT_NAME + '/Data/MCP_images/movie_scan_' + str(scan_num).zfill(3) + '.mp4'

        # subtract backgound
        if sub_bg == 1:
            imgs_data = self.img_data_subtract_bg(imgs_data)
            DIR_movie = './' + self.PROJECT_NAME + '/Data/MCP_images/movie_scan_' + str(scan_num).zfill(3) + '_subbg.mp4'

        # plotting
        print("Start generating images...")
        for snap_no in range(len(imgs_data)):
            fig = plt.figure(figsize=fig_size)
            plt.rcParams.update({'font.size': font_size})

            ax = fig.add_subplot(111)
            TITLE = 'snap ' + str(snap_no)
            if VARIABLE != '':
                TITLE = TITLE + ', ' + VARIABLE + ' = ' + str(round(var_data[snap_no], 2))
            ax.set_title(TITLE)
            img = plt.imshow(imgs_data[snap_no].T, origin='lower', vmin=v_min, vmax=v_max)
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

    def animation_mcp_imgs(self, imgs_data, v_range = [0, -1], fig_size = (10, 10), font_size = 20, show = 1, clean = 1):
        '''
            Generate animations with imgs_data.
            v_range = [v_min, v_max]
            v_max = -1: use highest count as vmax
        '''
        start_time = time.time()

        [v_min, v_max] = v_range
        if v_max == -1: # use highest count as vmax
            v_max = np.amax(imgs_data)

        # prepare directory
        DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/img_sub'
        if not os.path.exists(DIR):
            os.makedirs(DIR)
        DIR_movie = './' + self.PROJECT_NAME + '/Data/MCP_images/movie_imgs.mp4'

        # plotting
        print("Start generating images...")
        for snap_no in range(len(imgs_data)):
            fig = plt.figure(figsize=fig_size)
            plt.rcParams.update({'font.size': font_size})

            ax = fig.add_subplot(111)
            ax.set_title('snap ' + str(snap_no))
            img = plt.imshow(imgs_data[snap_no].T, origin='lower', vmin=v_min, vmax=v_max)
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

    def display_mcp_wireframe_1(self, img_data, frame_color = 'g', TITLE = '', z_max = -2, angle = (20, -60), fig_size = (5, 5), save = False, FILENAME = 'untitled'):
        '''
            Display mcp image by wireframe for img_data.
            z_max = -1: use highest count in imgs_data as z_max
            z_max = -2: use highest count in imgs_data[snap_shot] as z_max
        '''
        # plot
        fig = plt.figure(figsize=fig_size, dpi = 200)
        plt.rcParams.update({'font.size': 10})
        ax = fig.gca(projection='3d')

        # prepare data
        x = np.linspace(0, 1, self.MCP_pixel[0])
        y = np.linspace(0, 1, self.MCP_pixel[1])
        X, Y = np.meshgrid(x, y)
        Z = img_data.T

        # plot wireframe
        surf = ax.plot_wireframe(X, Y, Z, color=frame_color, linewidth=0.4, rcount = 64, ccount = 64)

        # set z_max
        if z_max == -2:
            z_max = np.max(img_data)
        ax.set_zlim(0, z_max)

        # clear ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

        # set axis color to white
        ax.w_xaxis.line.set_color('w')
        ax.w_yaxis.line.set_color('w')
        ax.w_zaxis.line.set_color('w')

        # set edge color to white
        ax.xaxis.pane.set_edgecolor('w')
        ax.yaxis.pane.set_edgecolor('w')
        ax.zaxis.pane.set_edgecolor('w')

        # clear grids
        ax.grid(False)

        # set rotation angle
        (Elev, Azim) = angle
        ax.view_init(elev=Elev, azim=Azim)

        if TITLE != '':
            plt.title(TITLE)
        plt.tight_layout()
        
        if save:
            # prepare directory
            DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/wireframe/'
            if not os.path.exists(DIR):
                os.makedirs(DIR)
            plt.savefig(DIR + 'wireframe_' + FILENAME + '.png', dpi = 200, format = 'png')

        plt.show()

        return

    def display_mcp_wireframe(self, imgs_data, snap_no, TITLE, z_max = -1, angle = (20, -60), fig_size = (5, 5), show = 1, save = 0, saveflag = 0):
        '''
            Display mcp image by wireframe for imgs_data[snap_no].
            z_max = -1: use highest count in imgs_data as z_max
            z_max = -2: use highest count in imgs_data[snap_shot] as z_max
        '''
        # get data
        img_data = imgs_data[snap_no]

        # plot
        fig = plt.figure(figsize=fig_size, dpi = 200)
        plt.rcParams.update({'font.size': 10})
        ax = fig.gca(projection='3d')

        # prepare data
        x = np.linspace(0, 1, self.MCP_pixel[0])
        y = np.linspace(0, 1, self.MCP_pixel[1])
        X, Y = np.meshgrid(x, y)
        Z = img_data.T

        # plot wireframe
        surf = ax.plot_wireframe(X, Y, Z, color='g', linewidth=0.4, rcount = 64, ccount = 64)

        # set z_max
        if z_max == -1:
            z_max = np.max(imgs_data)
        elif z_max == -2:
            z_max = np.max(img_data)
        ax.set_zlim(0, z_max)

        # clear ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

        # set axis color to white
        ax.w_xaxis.line.set_color('w')
        ax.w_yaxis.line.set_color('w')
        ax.w_zaxis.line.set_color('w')

        # set edge color to white
        ax.xaxis.pane.set_edgecolor('w')
        ax.yaxis.pane.set_edgecolor('w')
        ax.zaxis.pane.set_edgecolor('w')

        # clear grids
        ax.grid(False)

        # set rotation angle
        (Elev, Azim) = angle
        ax.view_init(elev=Elev, azim=Azim)

        plt.title(TITLE)
        plt.tight_layout()
        
        if save == 1:
            # prepare directory
            DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/wireframe/'
            if not os.path.exists(DIR):
                os.makedirs(DIR)
            plt.savefig(DIR + 'wireframe_' + str(saveflag).zfill(3) + '.png', dpi = 200, format = 'png')

        if show == 1:
            plt.show()
        elif show == 0:
            plt.close(fig)
        return

    def animation_mcp_wireframe(self, scan_num, TITLE_DETECTOR, TITLE_PREFIX, TITLE_SUFFIX, MOVIE_NAME = 'wireframe', order = 1, decimeter = 1, angle = (20, 30), fig_size = (5, 5), clean = 1, show = 1):
        '''
            Generate mcp animation with wireframe for scan No.scan_num.
            order =  1: positive sequence (as scan goes)
            order = -1: negative sequence

            2020.01.13
            to-do: using img_data instead of imgs_data when generating individual images.
        '''
        start_time = time.time()

        # get data
        imgs_data = self.img_scan_mcp(scan_num)
        variable = self.scan_value(scan_num, TITLE_DETECTOR)

        # prepare directory
        DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/wireframe'
        if not os.path.exists(DIR):
            os.makedirs(DIR)
        DIR_movie = './' + self.PROJECT_NAME + '/Data/MCP_images/movie_' + MOVIE_NAME + '.mp4'

        # generate plots
        print("Start generating images...")
        if order == 1:
            for i in range(len(imgs_data)):
                snap_no = i
                TITLE = TITLE_PREFIX + str(np.around(variable[snap_no], decimeter)) + TITLE_SUFFIX
                self.display_mcp_wireframe(imgs_data, snap_no, TITLE, angle = angle, fig_size = fig_size, show = 0, save = 1, saveflag = i)
        elif order == -1:
            for i in range(len(imgs_data)):
                snap_no = len(imgs_data) - i - 1
                TITLE = TITLE_PREFIX + str(np.around(variable[snap_no], decimeter)) + TITLE_SUFFIX
                self.display_mcp_wireframe(imgs_data, snap_no, TITLE, angle = angle, fig_size = fig_size, show = 0, save = 1, saveflag = i)
        print("Done.")

        # generate animation using ffmpeg
        print("Start generating animation...")
        os.system("ffmpeg -r 5 -i " + DIR + "/wireframe_%03d.png" 
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

    def animation_mcp_wireframe_sub(self, scan_nums, TITLE_DETECTOR, TITLE_PREFIX, TITLE_SUFFIX, MOVIE_NAME = 'wireframe', order = 1, decimeter = 1, angle = (20, 30), fig_size = (5, 5), clean = 1, show = 1):
        '''
            Generate mcp animation with wireframe for scan No.scan_nums[0] - No.scan_nums[1].
            order =  1: positive sequence (as scan goes)
            order = -1: negative sequence
        '''
        start_time = time.time()

        # get data
        imgs_data = self.img_scan_mcp(scan_nums[0]) - self.img_scan_mcp(scan_nums[1])
        variable = self.scan_value(scan_nums[0], TITLE_DETECTOR)

        # prepare directory
        DIR = './' + self.PROJECT_NAME + '/Data/MCP_images/wireframe'
        if not os.path.exists(DIR):
            os.makedirs(DIR)
        DIR_movie = './' + self.PROJECT_NAME + '/Data/MCP_images/movie_' + MOVIE_NAME + '.mp4'

        # generate plots
        print("Start generating images...")
        if order == 1:
            for i in range(len(imgs_data)):
                snap_no = i
                TITLE = TITLE_PREFIX + str(np.around(variable[snap_no], decimeter)) + TITLE_SUFFIX
                self.display_mcp_wireframe(imgs_data, snap_no, TITLE, angle = angle, fig_size = fig_size, show = 0, save = 1, saveflag = i)
        elif order == -1:
            for i in range(len(imgs_data)):
                snap_no = len(imgs_data) - i - 1
                TITLE = TITLE_PREFIX + str(np.around(variable[snap_no], decimeter)) + TITLE_SUFFIX
                self.display_mcp_wireframe(imgs_data, snap_no, TITLE, angle = angle, fig_size = fig_size, show = 0, save = 1, saveflag = i)
        print("Done.")

        # generate animation using ffmpeg
        print("Start generating animation...")
        os.system("ffmpeg -r 5 -i " + DIR + "/wireframe_%03d.png" 
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
            img = ax.scatter(db[0], db[1], db[2], c=db[3], cmap = colormap, linewidths = 0, s=marker_size)
        elif log == 1:
            img = ax.scatter(db[0], db[1], db[2], c=np.log(db[3]), cmap = colormap, linewidths = 0, s=marker_size)
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

    ################################################################################
    ############################## process SDD data ################################
    ################################################################################

    def data_scan_sdd(self, scan_num, print_info = 0):
        '''
            Read sdd data of scan_num from file, return energy & intensities data.
        '''
        FILENAME = self.scan_sdd_filename(scan_num)

        with open(FILENAME) as f:
            datablock = f.read()
            
        # extract basic info
        db = datablock.split("#C SDD Energy Scale\n")
        energy_bins_no = int(db[0].split('\n')[2].split(' ')[1])
        shutter_time = float(db[0].split('\n')[3].split(' ')[1])
        
        if print_info:
            print("========== Scan Info. ==========")
            print(db[0][0:-1])
            print("================================")

        # extract main data
        sdd_data_raw = db[1].split('\n#@MCA\n')
        energy_raw = sdd_data_raw[0].split('\n')[1:] # pop first line '#@CALIB -11.6 2.6351'
        intensities_raw = sdd_data_raw[1].split('\n')
        
        # map ASCII data to float data
        energy = np.array(list(map(float, energy_raw)))
        intensities = []
        for intensity_raw in intensities_raw:
            intensity = list(map(float, intensity_raw.split(' ')))
            intensities.append(intensity)
        intensities = np.transpose(intensities)

        # per second
        if self.ctps_key:
            intensities = intensities / shutter_time

        return energy, intensities

    def sdd_ROI(self, scan_num, ROI_range, print_info = 0):
        '''
            Return SDD ROI intensity with energy region between ROI_range[0] and ROI_range[1].
        '''
        energy, intensities = self.data_scan_sdd(scan_num)
        intensity = np.average(intensities[(energy >= ROI_range[0]) * (energy <= ROI_range[1])], axis = 0)
        return intensity

    ################################################################################
    ############################# visualize SDD data ###############################
    ################################################################################



    #######################################################################
    ############################## visualize ##############################
    #######################################################################

    def plot_scan(self, scan_nums, VARIABLE, DETECTOR, I0_BD3, linestyle, fig_size, font_size, TITLE = ''):
        plt.figure(figsize=fig_size)
        plt.rcParams.update({'font.size': font_size})
        for scan_num in scan_nums:
            x = self.scan_value(scan_num, VARIABLE)
            y = self.scan_value(scan_num, DETECTOR)
            if I0_BD3 == 1:
                y = y / self.scan_value(scan_num, 'I0_BD3')
            plt.plot(x, y, linestyle)
        legends = ['scan #' + legend for legend in list(map(str, scan_nums))]
        plt.xlabel(VARIABLE)
        plt.ylabel(DETECTOR)
        plt.legend(legends)
        plt.title(TITLE)
        plt.show()
        return

    def quick_plot(self, scan_nums, VARIABLE, DETECTOR, legends = '', smooth = 1, sort = 0, fig_size = (14, 7), font_size = 18, I0_BD3 = 1, xlabel = '', ylabel = '', TITLE = ''):
        '''
            A quick visualize of data.
        '''
        # prepare data
        x = []
        y = []
        for scan_num in scan_nums:
            x.append(self.scan_value(scan_num, VARIABLE))
            if I0_BD3 == 1:
                y.append(self.scan_value(scan_num, DETECTOR) / self.scan_value(scan_num, 'I0_BD3'))
            else:
                y.append(self.scan_value(scan_num, DETECTOR))

        # sort by variables
        if sort == 1:
            for i in range(len(y)):
                db = [x[i], y[i]]
                temp = np.transpose(db)
                db = np.transpose(temp[temp[:, 0].argsort()])
                [x[i], y[i]] = db

        x = np.array(x)
        y = np.array(y)

        # plot
        plt.figure(figsize=fig_size, dpi = 200)
        plt.rcParams.update({'font.size': font_size})
        colors = plt.cm.jet(np.linspace(0.1, 0.9, len(y)))

        if smooth == 0:
            for i in range(len(y)):
                plt.plot(x[i], y[i], 'o-', color = colors[i])
        else:
            for i in range(len(y)):
                plt.plot(x[i], y[i], 'o', color = colors[i])
            for i in range(len(y)):
                xx = np.linspace(x[i].min(), x[i].max(), 200)
                y_interp = interp1d(x[i], y[i], kind='cubic')
                y_smooth = savgol_filter(y_interp(xx), 13, 3)
                plt.plot(xx, y_smooth, '-', color = colors[i])

        if xlabel == '':
            xlabel = VARIABLE
        if ylabel == '':
            ylabel = DETECTOR
        
        #plt.colormaps
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if legends != '':
            plt.legend(legends)
        plt.title(TITLE)
        plt.show()
        return

    #############################################################################
    ################################### Trash ###################################
    #############################################################################

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
            raise ValueError('I0_BD3 must be 0 or 1!')

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
            raise ValueError('I0_BD3 must be 0 or 1!')

        plt.xlabel('tth')
        plt.ylabel(DETECTOR)
        plt.legend()
        plt.title(TITLE)
        plt.show()
        return
