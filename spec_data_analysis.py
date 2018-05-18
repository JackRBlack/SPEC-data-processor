def spec_scan_reader(FILEPATH, FILENAME, scan_num):
    
    '''
        Read SPEC data from a single data file with specific scan #. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            scan_num : [int] to specific the scan number.
            
        returns:
            [th, TEY, MCP, pm3] : [list]
            
        example:
            [th, TEY, MCP, pm3] = spec_scan_reader("./data/", "Sample_A", 56)
    ''' 
    
    import numpy as np
    
    PATH = FILEPATH + FILENAME + "_" + str(scan_num).zfill(3) + ".dat"
    data_origin = np.loadtxt(PATH)
    data = np.transpose(data_origin)
    th = data[0]
    TEY = data[7]
    MCP = data[8]
    pm3 = data[12]
    
    return [th, TEY, MCP, pm3]

def spec_Escan_reader(FILEPATH, FILENAME, scan_num):
    
    '''
        Read SPEC Escan data from a single data file with specific scan #. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            scan_num : [int] to specific the scan number.
            
        returns:
            [th, TEY, MCP, pm3] : [list]
            
        example:
            [th, TEY, MCP, pm3] = spec_scan_reader("./data/", "Sample_A", 56)
    ''' 
    
    import numpy as np
    
    PATH = FILEPATH + FILENAME + "_" + str(scan_num).zfill(3) + ".dat"
    data_origin = np.loadtxt(PATH)
    data = np.transpose(data_origin)
    E = data[0]
    TEY = data[6]
    MCP = data[7]
    pm3 = data[11]
    
    return [E, TEY, MCP, pm3]

def spec_data_reader(FILEPATH, FILENAME, scan_nums):
    
    '''
        Read SPEC data from mutiple data files with specific scan #s. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            scan_nums : [list] to specific the scan numbers.
            
        returns:
            [th, TEY, MCP, pm3] : [list]
            
        example:
            datablock = spec_data_reader("./data/", "Sample_A", [100, 101, 102])
    ''' 
    
    import numpy as np
    
    data = []
    
    for scan_num in scan_nums:
        data.append(spec_scan_reader(FILEPATH, FILENAME, scan_num))
    
    return data

def spec_plot(datablock, VARIABLE, DETECTOR, data_legends, LEGEND_PREFIX, LEGEND_SUFFIX):
    
    '''
        Plot a set of processed SPEC data together with legends. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            datablock : [list] data read from SPEC data file with function "spec_data_reader".
            DETECTOR : [string] must be chosen from "TEY", "MCP", "pm3".
            data_legends : [list] a set of legends, [20, 40, 60], e.g.
            LEGEND_PREFIX : [string] prefix to data_legends.
            LEGEND_SUFFIX : [string] suffix to data_legends.
            
        returns:
            a figure with multiple curves marked with legends
            
        example:
            spec_plot(datablock, "MCP", [100, 101, 102], "scan #", "")
            spec_plot(datablock, "MCP", [20, 50, 80], "T = ", " K")
    ''' 
    import matplotlib.pyplot as plt
    
    if DETECTOR == "TEY":
        det = 1
    elif DETECTOR == "MCP":
        det = 2
    elif DETECTOR == "pm3":
        det = 3
    else:
        raise ValueError('DETECTOR must be "TEY", "MCP" or "pm3"!')

    i = 0
    for legend in data_legends:
        plt.plot(datablock[i][0], datablock[i][det], label = LEGEND_PREFIX + str(legend) + LEGEND_SUFFIX)
        i = i + 1
    plt.xlabel(VARIABLE)
    plt.ylabel(DETECTOR)
    plt.legend()
    plt.show()
    return

def spec_data_average(FILEPATH, FILENAME, scan_nums):
    
    '''
        Read several SPEC scans' data and average them. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            scan_nums : [list] to specific the scan numbers.
            
        returns:
            [th, TEY, MCP, pm3] : [list]
            
        example:
            [th, TEY, MCP, pm3] = spec_data_average("./data/", "Sample_A", [100, 101, 102])
    ''' 
    
    import numpy as np
    
    data = []
    
    for scan_num in scan_nums:
        data.append(spec_scan_reader(FILEPATH, FILENAME, scan_num))
    
    th = data[0][0]
    TEY = data[0][1]
    MCP = data[0][2]
    pm3 = data[0][3]
    
    i = 1
    scan_nums_minus = scan_nums[:-1]
    for scan_num in scan_nums_minus:
        TEY = TEY + data[i][1]
        MCP = MCP + data[i][2]
        pm3 = pm3 + data[i][3]
        i = i + 1
        
    TEY = TEY / i
    MCP = MCP / i
    pm3 = pm3 / i
    
    return [th, TEY, MCP, pm3]

def spec_data_bg_subtracted_scans(FILEPATH, FILENAME, signal_scan_nums, bg_scan_nums):
    
    '''
        Subtract background from signals using SPEC data scan #s. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            signal_scan_nums : [list] to specific the scan numbers of signals.
            bg_scan_nums : [list] to specific the scan numbers of the background.
            
        returns:
            [th, TEY, MCP, pm3] : [list]
            
        example:
            [th, TEY, MCP, pm3] = spec_data_bg_subtracted_scans("./data/", "Sample_A", [100, 101, 102], [200, 201, 202])
    ''' 
    
    [th_sig, TEY_sig, MCP_sig, pm3_sig] = spec_data_average(FILEPATH, FILENAME, signal_scan_nums)
    [th_bg, TEY_bg, MCP_bg, pm3_bg] = spec_data_average(FILEPATH, FILENAME, bg_scan_nums)
    TEY = TEY_sig - TEY_bg
    MCP = MCP_sig - MCP_bg
    pm3 = pm3_sig - pm3_bg
    
    return [th_sig, TEY, MCP, pm3]

def spec_data_bg_subtracted(signal_data, bg_data):
    
    '''
        Subtract background from signals using readed/calculated SPEC data. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            signal_data : [list] = [th, TEY, MCP, pm3]
            bg_data : [list] = [th, TEY, MCP, pm3]
            
        returns:
            [th, TEY, MCP, pm3] : [list]
            
        example:
            [th, TEY, MCP, pm3] = spec_data_bg_subtracted(signal_data, bg_data)
    ''' 
    
    TEY = signal_data[1] - bg_data[1]
    MCP = signal_data[2] - bg_data[2]
    pm3 = signal_data[3] - bg_data[3]
    
    return [signal_data[0], TEY, MCP, pm3]