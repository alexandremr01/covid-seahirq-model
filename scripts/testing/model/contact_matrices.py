from .utils import read_matrix
import numpy as np
from scipy.sparse.linalg import eigs
from .nextgen import calculateR0

contact_matrix_all_file = 'contact_matrix_all.csv'
contact_matrix_home_file = 'contact_matrix_home.csv'
contact_matrix_work_file = 'contact_matrix_work.csv'
contact_matrix_school_file = 'contact_matrix_school.csv'
contact_matrix_other_file = 'contact_matrix_other.csv'

class ContactMatrixFactory:
    def __init__(self, age_strata):
        self.matrix_all = read_matrix(contact_matrix_all_file, age_strata)
        self.matrix_home = read_matrix(contact_matrix_home_file, age_strata)
        self.matrix_work = read_matrix(contact_matrix_work_file, age_strata)
        self.matrix_school = read_matrix(contact_matrix_school_file, age_strata)
        self.matrix_other = read_matrix(contact_matrix_other_file, age_strata)
        self.age_strata = age_strata


    def get_matrix(self, intervention, rel_pop, attack, age_strata, input_folder, parameters, R0):
        # Escala a matriz de contato pelo perfil demográfico e escalona a mesma pelo auto-valor dominante
        matrix_home = np.copy(self.matrix_home)
        matrix_work = np.copy(self.matrix_work)
        matrix_school = np.copy(self.matrix_school)
        matrix_other = np.copy(self.matrix_other)
        matrix = matrix_home + matrix_work + matrix_school + matrix_other

        attack_weight = [0.00253, 0.00268275, 0.00325998, 0.00587223, 0.00974128, 0.01562692, 0.02170502, 0.02817997,
                        0.03376366, 0.03846543, 0.04203209, 0.04679663, 0.05946109, 0.08191154, 0.0843, 0.0844]
        attack_weight = np.multiply(1.0, attack_weight)

        for i in range(0, age_strata):
            for j in range(0, age_strata):
                matrix[i, j] = attack_weight[i] * matrix[i, j] / rel_pop[i] #* pop[i] / pop[j]
                matrix_home[i, j] = attack_weight[i] * matrix_home[i, j] / rel_pop[i] #* pop[i] / pop[j]
                matrix_school[i, j] = attack_weight[i] * matrix_school[i, j] / rel_pop[i] #* pop[i] / pop[j]
                matrix_work[i, j] = attack_weight[i] * matrix_work[i, j] / rel_pop[i] #* pop[i] / pop[j]
                matrix_other[i, j] = attack_weight[i] * matrix_other[i, j] / rel_pop[i] #* pop[i] / pop[j]

        R0_raw = calculateR0(input_folder, matrix, parameters.xI[0], parameters.xA[0], False)
        R0_post = R0
        if attack is None:
            R0_ratio = R0/R0_raw
            R0_ratio_post = R0_post / R0_raw
        else:
            R0_ratio = attack
            R0_ratio_post = attack

        C_all_pre = R0_ratio * matrix
        C_home_pre = R0_ratio * matrix_home
        C_work_pre = R0_ratio * matrix_work
        C_school_pre = R0_ratio * matrix_school
        C_other_pre = R0_ratio * matrix_other

        C_home_post = R0_ratio_post * matrix_home
        C_work_post = R0_ratio_post * matrix_work
        C_school_post = R0_ratio_post * matrix_school
        C_other_post = R0_ratio_post * matrix_other
        C_all_post = C_home_post + C_work_post + C_school_post + C_other_post  ## itv_id = 1

        # Build matrix for scenarios
        I_old = np.diag(np.ones(age_strata))
        A_home = np.diag(np.ones(age_strata))
        B_other = np.diag(np.ones(age_strata))
        B_school = np.diag(np.ones(age_strata))
        B_lock = np.diag(np.ones(age_strata))
        B_strong_lock = np.diag(np.ones(age_strata))
        W_work = np.diag(np.ones(age_strata))
        W_lock = np.diag(np.ones(age_strata))
        W_strong_lock = np.diag(np.ones(age_strata))
        for i in range(age_strata - 5, age_strata):
            I_old[i, i] = 0.1
        for i in range(0, 4):
            A_home[i, i] = 1.5
        for i in range(4, age_strata):
            A_home[i, i] = 1.1
        for i in range(0, 4):
            B_school[i, i] = 0.4
        for i in range(0, 4):
            B_other[i, i] = 0.4
        for i in range(4, age_strata):
            B_other[i, i] = 0.6
        for i in range(0, 4):
            B_lock[i, i] = 0.3
        for i in range(4, age_strata):
            B_lock[i, i] = 0.5
        for i in range(0, 4):
            B_strong_lock[i, i] = 0.1
        for i in range(4, age_strata):
            B_strong_lock[i, i] = 0.1
        for i in range(0, age_strata):
            W_work[i, i] = 0.6
        for i in range(0, age_strata):
            W_lock[i, i] = 0.4
        for i in range(0, age_strata):
            W_strong_lock[i, i] = 0.3

        C_work_old = np.dot(np.dot(I_old, C_work_post), I_old)
        C_school_old = np.dot(np.dot(I_old, C_school_post), I_old)
        C_other_old = np.dot(np.dot(I_old, C_other_post), I_old)

        if intervention == 0:
            return C_all_pre
        elif intervention == 1:
            return C_all_post
        elif intervention == 2:
            # Fechamento de escolas apenas
            C_all_school = np.dot(A_home, C_home_post) + C_work_post + np.dot(B_school, C_other_post)
            return C_all_school
        elif intervention == 3:
            # Fechamento de escolas e distanciamento social
            C_all_school_other = np.dot(A_home, C_home_post) + C_work_post + + np.dot(B_other, C_other_post)
            return C_all_school_other
        elif intervention == 4:
            # Fechamento de escolas, distanciamento social e trabalho
            C_all_school_other_work = np.dot(A_home, C_home_post) + np.dot(W_work, C_work_post) + np.dot(B_other, C_other_post)
            return C_all_school_other_work
        elif intervention == 5:
            # Lockdown sem isolamento de idoso
            C_all_lock = np.dot(A_home, C_home_post) + np.dot(W_lock, C_work_post) + np.dot(B_lock, C_other_post)
            return C_all_lock
        elif intervention == 6:
            # Isolamento dos idosos em, com redução de X% dos seus contatos
            C_all_old = C_home_post + C_work_old + C_school_old + C_other_old
            return C_all_old
        elif intervention == 7:
            # Isolamento dos idosos, fechamento de escolas
            C_all_old_school = np.dot(A_home, C_home_post) + C_work_old + np.dot(B_school, C_other_old)
            return C_all_old_school
        elif intervention == 8:
            # Isolamento dos idosos, fechamento de escolas e distanciamento social
            C_all_old_school_other = np.dot(A_home, C_home_post) + C_work_old + np.dot(B_other, C_other_old)
            return C_all_old_school_other
        elif intervention == 9:
            # Isolamento dos idosos, fechamento de escolas, distanciamento social e distanciamento no trabalho
            C_all_old_school_other_work = np.dot(A_home, C_home_post) + np.dot(W_work, C_work_old) + np.dot(B_other, C_other_old)
            return C_all_old_school_other_work
        elif intervention == 10:
            # Lockdown com isolamento de idoso
            C_all_old_lock = np.dot(A_home, C_home_post) + np.dot(W_lock, C_work_old) + np.dot(B_lock, C_other_old)
            return C_all_old_lock
        elif intervention == 11:
            # Lockdown com isolamento de idoso
            C_all_old_strong_lock = np.dot(A_home, C_home_post) + np.dot(W_strong_lock, C_work_old) + np.dot(B_strong_lock,
                                                                                                        C_other_old)
            return C_all_old_strong_lock
        elif intervention == 12:
             # Apenas distanciamento social
            C_all_other = C_home_post + C_work_post + C_school_post + np.dot(B_other, C_other_post)
            return C_all_other
        elif intervention == 13:
            # Isolamento de idoso + distanciamento social
            C_all_old_other = C_home_post + C_work_old + C_school_old + np.dot(B_other, C_other_old)
            return C_all_old_other