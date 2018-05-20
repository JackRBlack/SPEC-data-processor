# read data

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
            [th, I0_BD3, TEY, MCP, pm3] : [list]
            
        example:
            [th, I0_BD3, TEY, MCP, pm3] = spec_scan_reader("./data/", "Sample_A", 56)
    ''' 
    
    import numpy as np
    
    PATH = FILEPATH + FILENAME + "_" + str(scan_num).zfill(3) + ".dat"
    data_origin = np.loadtxt(PATH)
    data = np.transpose(data_origin)
    th = data[0]
    I0_BD3 = data[6]
    TEY = data[7]
    MCP = data[8]
    pm3 = data[12]
    
    return [th, I0_BD3, TEY, MCP, pm3]

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
            [E, I0_BD3, TEY, MCP, pm3] : [list]
            
        example:
            [E, I0_BD3, TEY, MCP, pm3] = spec_Escan_reader("./data/", "Sample_A", 56)
    ''' 
    
    import numpy as np
    
    PATH = FILEPATH + FILENAME + "_" + str(scan_num).zfill(3) + ".dat"
    data_origin = np.loadtxt(PATH)
    data = np.transpose(data_origin)
    E = data[0]
    I0_BD3 = data[5]
    TEY = data[6]
    MCP = data[7]
    pm3 = data[11]
    
    return [E, I0_BD3, TEY, MCP, pm3]

def spec_th2th_reader(FILEPATH, FILENAME, scan_num):
    
    '''
        Read SPEC Escan data from a single data file with specific scan #. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            scan_num : [int] to specific the scan number.
            
        returns:
            [tth, th, I0_BD3, TEY, MCP, pm3] : [list]
            
        example:
            [tth, th, I0_BD3, TEY, MCP, pm3] = spec_th2th_reader("./data/", "Sample_A", 56)
    ''' 
    
    import numpy as np
    
    PATH = FILEPATH + FILENAME + "_" + str(scan_num).zfill(3) + ".dat"
    data_origin = np.loadtxt(PATH)
    data = np.transpose(data_origin)
    tth = data[0]
    th = data[1]
    I0_BD3 = data[7]
    TEY = data[8]
    MCP = data[9]
    pm3 = data[13]
    
    return [tth, th, I0_BD3, TEY, MCP, pm3]

def spec_data_reader(SCAN_FORMAT, FILEPATH, FILENAME, scan_nums):
    
    '''
        Read SPEC data from mutiple data files with specific scan #s. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            SCAN_FORMAT : [string] to indicate the scan format, "ascan" or "a2scan".
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            scan_nums : [list] to specific the scan numbers.
            
        returns:
            [th, I0_BD3, TEY, MCP, pm3] or [tth, th, I0_BD3, TEY, MCP, pm3] : [list]
            
        example:
            datablock = spec_data_reader("./data/", "Sample_A", [100, 101, 102])
    ''' 
    
    import numpy as np
    
    data = []
    
    if SCAN_FORMAT == "ascan":
        for scan_num in scan_nums:
            data.append(spec_scan_reader(FILEPATH, FILENAME, scan_num))
    elif SCAN_FORMAT == "a2scan":
        for scan_num in scan_nums:
            data.append(spec_th2th_reader(FILEPATH, FILENAME, scan_num))
    else:
        raise ValueError('SCAN_FORMAT must be "ascan" or "a2scan"!')
    return data

# data process

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
            [th, IO_BD3, TEY, MCP, pm3] : [list]
            
        example:
            [th, IO_BD3, TEY, MCP, pm3] = spec_data_average("./data/", "Sample_A", [100, 101, 102])
    ''' 
    
    import numpy as np
    
    data = []
    
    for scan_num in scan_nums:
        data.append(spec_scan_reader(FILEPATH, FILENAME, scan_num))
    
    th = data[0][0]
    I0_BD3 = data[0][1]
    TEY = data[0][2] / data[0][1]
    MCP = data[0][3] / data[0][1]
    pm3 = data[0][4] / data[0][1]
    
    i = 1
    scan_nums_minus = scan_nums[:-1]
    for scan_num in scan_nums_minus:
        TEY = TEY + data[i][2] / data[i][1]
        MCP = MCP + data[i][3] / data[i][1]
        pm3 = pm3 + data[i][4] / data[i][1]
        i = i + 1
        
    TEY = TEY / i
    MCP = MCP / i
    pm3 = pm3 / i
    
    return [th, I0_BD3, TEY, MCP, pm3]

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
    
    [th_sig, I0_BD3_sig, TEY_sig, MCP_sig, pm3_sig] = spec_data_average(FILEPATH, FILENAME, signal_scan_nums)
    [th_bg, I0_BD3_bg, TEY_bg, MCP_bg, pm3_bg] = spec_data_average(FILEPATH, FILENAME, bg_scan_nums)
    TEY = TEY_sig - TEY_bg
    MCP = MCP_sig - MCP_bg
    pm3 = pm3_sig - pm3_bg
    
    return [th_sig, I0_BD3_sig, TEY, MCP, pm3]

def spec_data_bg_subtracted(signal_data, bg_data):
    
    '''
        Subtract background from signals using readed/calculated SPEC data. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            signal_data : [list] = [th, TEY, MCP, pm3]
            bg_data : [list] = [th, TEY, MCP, pm3]
            
        returns:
            [th, I0_BD3, TEY, MCP, pm3] : [list]
            
        example:
            [th, I0_BD3, TEY, MCP, pm3] = spec_data_bg_subtracted(signal_data, bg_data)
    ''' 
    
    TEY = signal_data[2] - bg_data[2]
    MCP = signal_data[3] - bg_data[3]
    pm3 = signal_data[4] - bg_data[4]
    
    return [signal_data[0], signal_data[1], TEY, MCP, pm3]

# generate figure

def spec_plot(datablock, SCAN_FORMAT, VARIABLE, DETECTOR, I0_BD3, data_legends, LEGEND_PREFIX, LEGEND_SUFFIX, TITLE):
    
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
    import matplotlib.pyplot as plt
    
    plt.figure()
    
    if SCAN_FORMAT == "ascan":
        
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
    
    elif SCAN_FORMAT == "a2scan":
        
        if DETECTOR == "TEY":
            det = 3
        elif DETECTOR == "MCP":
            det = 4
        elif DETECTOR == "pm3":
            det = 5
        else:
            raise ValueError('DETECTOR must be "TEY", "MCP" or "pm3"!')
        
        if VARIABLE == "th":
            if I0_BD3 == 0:
                i = 0
                for legend in data_legends:
                    plt.plot(datablock[i][1], datablock[i][det], label = LEGEND_PREFIX + str(legend) + LEGEND_SUFFIX)
                    i = i + 1
            elif I0_BD3 == 1:
                i = 0
                for legend in data_legends:
                    plt.plot(datablock[i][1], datablock[i][det]/datablock[i][2], label = LEGEND_PREFIX + str(legend) + LEGEND_SUFFIX)
                    i = i + 1
            else:
                raise ValueError('IO_BD3 must be 0 or 1!')
        elif VARIABLE == "tth":
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
        else:
            raise ValueError('In a2scan VARIABLE must be "th" or "tth"!')
    plt.xlabel(VARIABLE)
    plt.ylabel(DETECTOR)
    plt.legend()
    plt.title(TITLE)
    plt.show()
    return

def spec_scatter3_th2th(datablock, DETECTOR, I0_BD3, scan_nums, TITLE):
    '''
        Scatter 3D-data: (th, tth, counts). 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            datablock : [[tth, th, I0_BD3, TEY, MCP, pm3], ...] data read from SPEC data file with function "spec_data_reader".
            DETECTOR : [string] must be chosen from "TEY", "MCP", "pm3".
            I0_BD3 : [int] 0 - do not divided by I0_BD3, 1 - divided by I0_BD3.
            scan_nums : [list] to specific the scan numbers.
            TITLE : [string] figure title
            
        returns:
            a 3D figure with lots of dots
            
        example:
            spec_plot(datablock, "a2scan", "tth", "MCP", 1, [100, 101, 102], "scan #", "")
            spec_plot(datablock, "ascan", "th", "MCP", 0, [20, 50, 80], "T = ", " K")
    ''' 
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    if DETECTOR == "TEY":
        det = 3
    elif DETECTOR == "MCP":
        det = 4
    elif DETECTOR == "pm3":
        det = 5
    else:
        raise ValueError('DETECTOR must be "TEY", "MCP" or "pm3"!')
            
    if I0_BD3 == 0:
        i = 0
        for scan_num in scan_nums:
            ax.scatter(datablock[i][0], datablock[i][1], datablock[i][det], s=20, c=None, depthshade=True)
            i = i + 1
    elif I0_BD3 == 1:
        i = 0
        for scan_num in scan_nums:
            ax.scatter(datablock[i][0], datablock[i][1], datablock[i][det]/datablock[i][2], s=20, c=None, depthshade=True)
            i = i + 1
    else:
        raise ValueError('IO_BD3 must be 0 or 1!')
    ax.set_xlabel('tth')
    ax.set_ylabel('th')
    ax.set_zlabel(DETECTOR)
    plt.title(TITLE)
    plt.show()

def spec_show_a2scan_route(FILEPATH, FILENAME, scan_nums, line_color):
    '''
        Draw th2th scan route. 
        
        Author: Wenjie Chen 
        E-mail: wenjiechen@pku.edu.cn
        
        args:
            FILEPATH : [string] to indicate where the data file is located.
            FILENAME : [string] to specific the data file name.
            scan_nums : [list] to specific the scan numbers.
            line_color : [tupple] to define the line color (RGB), for example (0.1,0.5,0.1).
            
        returns:
            Serveral lines in current canvas.
            
        example:
            plt.figure()
            spec_show_a2scan_route("./data/", "sample_a", [123, 124, 125], (0.1,0.5,0.1))
            plt.show()
    '''
    
    for scan_num in scan_nums:
        data = spec_th2th_reader(FILEPATH, FILENAME, scan_num)
        tth = [data[0][0], data[0][-1]]
        th = [data[1][0], data[1][-1]]
        plt.plot(th, tth, c=line_color)