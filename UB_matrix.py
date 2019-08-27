########## INFO ##########
print("########################################")
print("Project: UB_matrix")
print("Version: 1.0")
print("Last Update: 2019.08.27")
print("----------------------------------------")
print("Author: Wenjie Chen")
print("E-mail: wenjiechen@pku.edu.cn")
print("########################################")
##########################

import numpy as np

class crystal:
    def __init__(self):
        self.lattice_lengths = np.zeros([3]) # in Angstrom
        self.lattice_angles = np.zeros([3]) # in degrees
        self.reciprocal_lattice_lengths = np.zeros([3]) # in Angstrom^-1, physicist notations (with 2pi)
        self.reciprocal_lattice_angles = np.zeros([3]) # in degrees

    def info(self):
        print("========== Lattice Constant ==========")
        print("Lattice:")
        print(f"[a, b, c] = {self.lattice_constant} Angstrom")
        print(f"[alpha, beta, gamma] = {self.lattice_angles} Degrees")
        print("--------------------------------------")
        print("Reciprocal lattice:")
        print(f"[a*, b*, c*] = {self.reciprocal_lattice_constant} Angstrom^-1")
        print(f"[alpha*, beta*, gamma*] = {self.reciprocal_lattice_angles} Degrees")
        return

    def set_lat(self, lengths, angles):
        '''
            Set lattice constant and calculate reciprocal lattice constant.
        '''
        # set lattice constant
        self.lattice_lengths = lengths
        self.lattice_angles = angles

        # calculate reciprocal lattice constant
        thetas = np.deg2rad(self.lattice_angles) # convert degrees to rads
        V_angles = np.sqrt(1 - np.cos(thetas[0])**2 - np.cos(thetas[1])**2 - np.cos(thetas[2])**2 + \
                   2 * np.cos(thetas[0]) * np.cos(thetas[1]) * np.cos(thetas[2])) # V = V_angles * a * b * c

        self.reciprocal_lattice_lengths[0] = np.sin(thetas[0]) / self.lattice_lengths[0] / V_angles * 2 * np.pi
        self.reciprocal_lattice_lengths[1] = np.sin(thetas[1]) / self.lattice_lengths[1] / V_angles * 2 * np.pi
        self.reciprocal_lattice_lengths[2] = np.sin(thetas[2]) / self.lattice_lengths[2] / V_angles * 2 * np.pi

        self.reciprocal_lattice_angles[0] = np.arccos((np.cos(thetas[1]) * np.cos(thetas[2]) - np.cos(thetas[0])) / (np.sin(thetas[1]) * np.sin(thetas[2])))
        self.reciprocal_lattice_angles[1] = np.arccos((np.cos(thetas[0]) * np.cos(thetas[2]) - np.cos(thetas[1])) / (np.sin(thetas[0]) * np.sin(thetas[2])))
        self.reciprocal_lattice_angles[2] = np.arccos((np.cos(thetas[0]) * np.cos(thetas[1]) - np.cos(thetas[2])) / (np.sin(thetas[0]) * np.sin(thetas[1])))
        self.reciprocal_lattice_angles = np.rad2deg(self.reciprocal_lattice_angles)
        return

class beam:
    def __init__(self):
        self.beam_energy = 0 # in eV
        self.beam_wavelength = 0 # in Angstrom
        self.beam_k = 0 # in Angstrom^-1

    def info(self):
        print("=========== Beam Constant ============")
        print(f"energy = {self.beam_energy} eV")
        print(f"wavelength = {self.beam_wavelength} Angstrom")
        print(f"k = {self.beam_k} Angstrom^-1")
        return
    
    def eV2k(self, energy):
        k = 2 * np.pi / (12398.42 / energy)
        return k

    def set_beam(self, energy):
        '''
            Set beam energy and calculate relavant constant.
        '''
        # set beam energy
        self.beam_energy = energy

        # calculate other beam constant
        self.beam_wavelength = 12398.42 / self.beam_energy
        self.beam_k = 2 * np.pi / self.beam_wavelength
        return

class scatter:
    '''
        Functions:

            Q ----B_matrix----> Q_c ----U_matrix----> Q_l

            - Q = [h, k, l]

            - Q_c = Q in the cartesian coordinate system tied to the sample

            - Q_l = k_i - k_f which is Q in the laboratory frame

        Conventions:

            - incident beam is along z axis, with x in the horizontal plane, and y pointing up.
            - the scattered beam makes an angle theta with the z axis, and its projection on the 
              (x, y) plane makes an angle phi with the x axis.

                Q_l = [- k_f sin(theta) cos(phi), - k_f sin(theta) sin(phi), k_i - k_f cos(theta)]
    '''

    def __init__(self):
        self.crystal = crystal()
        self.beam = beam()

        self.peak_1 = np.zeros([3, 3]) # [[h, k, l], [theta, phi, energy], [omega, mu, nu]]
        self.peak_2 = np.zeros([3, 3]) # [[h, k, l], [theta, phi, energy], [omega, mu, nu]]

        self.B_matrix = np.zeros([3, 3])
        self.U_matrix = np.zeros([3, 3])
        self.UB_matrix = np.zeros([3, 3])
        self.UB_matrix_inv = np.zeros([3, 3])

    def info(self):
        self.crystal.info()
        self.beam.info()
        return

    def set_lat(self, lengths, angles):
        self.crystal.set_lat(lengths, angles)
        return

    def set_beam(self, energy):
        self.beam.set_beam(energy)
        return

    def set_peak_1(self, index, angles, rotations):
        self.peak_1[0] = index
        self.peak_1[1] = angles
        self.peak_1[2] = rotations
        return

    def set_peak_2(self, index, angles, rotations):
        self.peak_2[0] = index
        self.peak_2[1] = angles
        self.peak_2[2] = rotations
        return

    def angles2Q(self, theta, phi, energy):
        '''
            Convert detector position angles to Q_l in the laboratory frame.
            theta: angle between k_f and k_i
            phi: angle between the projection of k_f on the (x, y) plane and the x axis
            angles in degrees.
        '''
        k_i = self.beam.beam_k
        k_f = self.beam.eV2k(energy)

        Q_l = np.zeros([3])
        Q_l[0] = - k_f * np.sin(np.deg2rad(theta)) * np.cos(np.deg2rad(phi))
        Q_l[1] = - k_f * np.sin(np.deg2rad(theta)) * np.sin(np.deg2rad(phi))
        Q_l[2] = k_i - k_f * np.cos(np.deg2rad(theta))
        return Q_l

    def R_1(self, omega):
        '''
            Return rotation matrix which rotates along y axis.
            omega in degrees.
        '''
        M = np.zeros([3, 3])
        M[1][1] = 1
        M[0][0] = np.cos(np.deg2rad(omega))
        M[0][2] = np.sin(np.deg2rad(omega))
        M[2][0] = - np.sin(np.deg2rad(omega))
        M[2][2] = np.cos(np.deg2rad(omega))
        return M

    def R_2(self, mu):
        '''
            Return rotation matrix which rotates along z axis.
            mu in degrees.
        '''
        M = np.zeros([3, 3])
        M[2][2] = 1
        M[0][0] = np.cos(np.deg2rad(mu))
        M[0][1] = - np.sin(np.deg2rad(mu))
        M[1][0] = np.sin(np.deg2rad(mu))
        M[1][1] = np.cos(np.deg2rad(mu))
        return M

    def R_3(self, nu):
        '''
            Return rotation matrix which rotates along x axis.
            mu in degrees.
        '''
        M = np.zeros([3, 3])
        M[0][0] = 1
        M[1][1] = np.cos(np.deg2rad(nu))
        M[1][2] = - np.sin(np.deg2rad(nu))
        M[2][1] = np.sin(np.deg2rad(nu))
        M[2][2] = np.cos(np.deg2rad(nu))
        return M

    def R(self, omega, mu, nu):
        '''
            Return rotation matrix which rotates along 3 axes.
            angles in degrees.
        '''
        M = np.dot(self.R_1(omega), np.dot(self.R_2(mu), self.R_3(nu)))
        return M

    def cal_B_matrix(self):
        '''
            Calculate the B matrix,
            which converts (h, k, l) to sample orthogonal reference frame.

            Project the following vector

                Q_c / 2 pi = h a* + k b* + l c*

            to the sample orthogonal reference frame, where

                i = a* unit vector
                j = perpendicular to i, in the (a* b*) plane
                k = perpendicular to i * j
        '''
        self.B_matrix[0][0] =  self.crystal.reciprocal_lattice_lengths[0]
        self.B_matrix[0][1] =  self.crystal.reciprocal_lattice_lengths[1] * np.cos(np.deg2rad(self.crystal.reciprocal_lattice_angles[2]))
        self.B_matrix[0][2] =  self.crystal.reciprocal_lattice_lengths[2] * np.cos(np.deg2rad(self.crystal.reciprocal_lattice_angles[1]))

        self.B_matrix[1][1] =  self.crystal.reciprocal_lattice_lengths[1] * np.sin(np.deg2rad(self.crystal.reciprocal_lattice_angles[2]))
        self.B_matrix[1][2] =  self.crystal.reciprocal_lattice_lengths[2] * np.sin(np.deg2rad(self.crystal.reciprocal_lattice_angles[1])) * \
                               - np.cos(np.deg2rad(self.crystal.lattice_angles[0]))

        self.B_matrix[2][2] =  2 * np.pi / self.crystal.lattice_lengths[2] # no 2pi in crystallographic notations 
        return

    def cal_U_matrix(self):
        '''
            Calculate the U matrix using two Bragg peaks.
        '''
        # prepare peak 1
        v_1 = self.peak_1[0]
        Q_1 = self.angles2Q(self.peak_1[1][0], self.peak_1[1][1], self.peak_1[1][2])
        R_1 = self.R(self.peak_1[2][0], self.peak_1[2][1], self.peak_1[2][2])
        Bv_1 = np.dot(self.B_matrix, v_1)
        RQ_1 = np.dot(np.transpose(R_1), Q_1) # transpose(R) = inv(R)

        # prepare peak 2
        v_2 = self.peak_2[0]
        Q_2 = self.angles2Q(self.peak_2[1][0], self.peak_2[1][1], self.peak_2[1][2])
        R_2 = self.R(self.peak_2[2][0], self.peak_2[2][1], self.peak_2[2][2])
        Bv_2 = np.dot(self.B_matrix, v_2)
        RQ_2 = np.dot(np.transpose(R_2), Q_2) # transpose(R) = inv(R)

        # calculate T_c matrix
        T_c = np.zeros([3, 3])

        T_c[0] = Bv_1 / np.sqrt(np.dot(Bv_1, Bv_1))
        temp = np.cross(Bv_1, Bv_2)
        T_c[2] = temp / np.sqrt(np.dot(temp, temp))
        T_c[1] = np.cross(T_c[2], T_c[0])
        T_c = np.transpose(T_c)

        # calculate T_nu matrix
        T_nu = np.zeros([3, 3])

        T_nu[0] = RQ_1 / np.sqrt(np.dot(RQ_1, RQ_1))
        temp = np.cross(RQ_1, RQ_2)
        T_nu[2] = temp / np.sqrt(np.dot(temp, temp))
        T_nu[1] = np.cross(T_nu[2], T_nu[0])
        T_nu = np.transpose(T_nu)

        # calculate U matrix
        self.U_matrix = np.dot(T_nu, np.transpose(T_c)) # transpose(T_c) = inv(T_c)
        
        return

    def cal_UB_matrix(self):
        '''
            Calculate UB matrix.
        '''
        self.cal_B_matrix()
        self.cal_U_matrix()
        self.UB_matrix = np.dot(self.U_matrix, self.B_matrix)
        self.UB_matrix_inv = np.linalg.inv(self.UB_matrix)
        return

    def position2hkl(self, angles, rotations):
        Q = self.angles2Q(angles[0], angles[1], angles[2])
        R = self.R(rotations[0], rotations[1], rotations[2])

        hkl = np.dot(self.UB_matrix_inv, np.dot(np.transpose(R), Q))

        return hkl

    def hkl2Q(self, hkl, rotations):
        R = self.R(rotations[0], rotations[1], rotations[2])
        Q = np.dot(R, np.dot(self.UB_matrix, hkl))
        return Q