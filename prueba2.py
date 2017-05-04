#!/usr/bin/env python
# coding: utf-8


class Structures:
    """
    Resolves structures by the stiffness method

    """

    def __init__(self, joints, pole, internal_loads, external_loads):
        self.joints = joints
        self.pole = pole
        self.internal_loads = internal_loads
        self.external_loads = external_loads

    @staticmethod
    def zeros(height, width):
        zero = []
        for i in range(height):
            zero1 = []
            for j in range(width):
                zero1.append(0.0)
            zero.append(zero1)
        return zero

    @staticmethod
    def transpose(m):
        m1 = Structures.zeros(len(m), len(m[0]))
        for i in range(len(m)):
            for j in range(len(m[0])):
                m1[j][i] = m[i][j]
        return m1

    @staticmethod
    def matrix_solve(m, c):
        """For L. Code by: Michael Halls-Moore"""
        n = len(m)
        l = Structures.zeros(n, n)
        for i in range(n):
            for k in range(i + 1):
                tmp_sum = sum(l[i][j] * l[k][j] for j in range(k))

                if i == k:
                    l[i][k] = (m[i][i] - tmp_sum) ** 0.5
                else:
                    l[i][k] = (1.0 / l[k][k] * (m[i][k] - tmp_sum))
        z = Structures.zeros(n, 1)
        for i in range(n):
            sum1 = 0
            for k in range(i):
                sum1 = sum1 + l[i][k] * z[k][0]
            z[i][0] = (c[i][0] - sum1) / l[i][i]
        x = Structures.zeros(n, 1)
        for i in range(n - 1, -1, -1):
            sum1 = 0
            for k in range(i + 1, n):
                sum1 = sum1 + l[k][i] * x[k][0]
            x[i][0] = (z[i][0] - sum1) / l[i][i]
        return x

    @staticmethod
    def scalar_product(m, n):
        for i in range(len(m)):
            for j in range(len(m[0])):
                m[i][j] = m[i][j] * n
        return m

    @staticmethod
    def remove_column(m, n):
        for i in range(len(m)):
            m[i].pop(n - 1)
        m1 = Structures.scalar_product(m, 1.0)
        return m1

    @staticmethod
    def remove_row(m, n):
        m.pop(n - 1)
        m1 = Structures.scalar_product(m, 1.0)
        return m1

    @staticmethod
    def sum_matrix(m1, m2):
        m3 = Structures.zeros(len(m1), len(m1[0]))
        for row in range(len(m1)):
            for column in range(len(m1[0])):
                m3[row][column] = m1[row][column] + m2[row][column]
        return m3

    @staticmethod
    def longitude(joints, pole):
        l = []
        for bars in range(len(pole)):
            l.append(((joints[pole[bars][1] - 1][0] - joints[pole[bars][0] - 1][0]) ** 2 +
                      (joints[pole[bars][1] - 1][1] - joints[pole[bars][0] - 1][1]) ** 2) ** 0.5)
        return l

    @staticmethod
    def moment_inertia(pole):
        inertia = []
        for bars in range(len(pole)):
            inertia.append(pole[bars][2]*pole[bars][3]**3/12)
        return inertia

    @staticmethod
    def area_poles(pole):
        area = []
        for bars in range(len(pole)):
            area.append(pole[bars][2] * pole[bars][3])
        return area

    @staticmethod
    def global_stiffness(joints, pole):
        stiff = []
        n_joints = len(joints)
        longs = Structures.longitude(joints, pole)
        model = [bars[4] for bars in pole]
        inertia = Structures.moment_inertia(pole)
        area = Structures.area_poles(pole)
        for i in range(len(pole)):
            cos_x = (joints[pole[i][1] - 1][0] - joints[pole[i][0] - 1][0]) / longs[i]
            cos_y = (joints[pole[i][1] - 1][1] - joints[pole[i][0] - 1][1]) / longs[i]
            k = Structures.zeros(3*n_joints, 3*n_joints)
            k[pole[i][0] * 3 - 3][pole[i][0] * 3 - 3] = (cos_x ** 2 * model[i] * area[i] / longs[
                i] + cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][0] * 3 - 3][pole[i][0] * 3 - 2] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][0] * 3 - 2][pole[i][0] * 3 - 3] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][0] * 3 - 2][pole[i][0] * 3 - 2] = (cos_y ** 2 * model[i] * area[i] / longs[
                i] + cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][0] * 3 - 3][pole[i][0] * 3 - 1] = (-6 * cos_y * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 1][pole[i][0] * 3 - 3] = (-6 * cos_y * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 1][pole[i][0] * 3 - 1] = (4 * model[i] * inertia[i] / longs[i])
            k[pole[i][0] * 3 - 2][pole[i][0] * 3 - 1] = (6 * cos_x * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 1][pole[i][0] * 3 - 2] = (6 * cos_x * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 3][pole[i][1] * 3 - 3] = - (cos_x ** 2 * model[i] * area[i] / longs[
                i] + cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][0] * 3 - 3][pole[i][1] * 3 - 2] = - ((model[i] * area[i] / longs[i] - 12 * model[
                i] * inertia[i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][0] * 3 - 2][pole[i][1] * 3 - 3] = - ((model[i] * area[i] / longs[i] - 12 * model[
                i] * inertia[i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][0] * 3 - 2][pole[i][1] * 3 - 2] = - (cos_y ** 2 * model[i] * area[i] / longs[
                i] + cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][0] * 3 - 3][pole[i][1] * 3 - 1] = (-6 * cos_y * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 1][pole[i][1] * 3 - 3] = (6 * cos_y * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 2][pole[i][1] * 3 - 1] = (6 * cos_x * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 1][pole[i][1] * 3 - 2] = (-6 * cos_x * model[i] * inertia[i] / longs[
                i] ** 2)
            k[pole[i][0] * 3 - 1][pole[i][1] * 3 - 1] = (2 * model[i] * inertia[i] / longs[i])
            k[pole[i][1] * 3 - 3][pole[i][0] * 3 - 3] = - (cos_x ** 2 * model[i] * area[i] / longs[
                i] + cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][1] * 3 - 2][pole[i][0] * 3 - 3] = - ((model[i] * area[i] / longs[i] - 12 * model[
                i] * inertia[i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][1] * 3 - 3][pole[i][0] * 3 - 2] = - ((model[i] * area[i] / longs[i] - 12 * model[
                i] * inertia[i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][1] * 3 - 2][pole[i][0] * 3 - 2] = - (cos_y ** 2 * model[i] * area[i] / longs[
                i] + cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][1] * 3 - 1][pole[i][0] * 3 - 3] = (-6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[pole[i][1] * 3 - 3][pole[i][0] * 3 - 1] = (6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[pole[i][1] * 3 - 1][pole[i][0] * 3 - 2] = (6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[pole[i][1] * 3 - 2][pole[i][0] * 3 - 1] = (-6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[pole[i][1] * 3 - 1][pole[i][0] * 3 - 1] = (2 * model[i] * inertia[i] / longs[i])
            k[pole[i][1] * 3 - 3][pole[i][1] * 3 - 3] = (cos_x ** 2 * model[i] * area[i] / longs[
                i] + cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][1] * 3 - 3][pole[i][1] * 3 - 2] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][1] * 3 - 2][pole[i][1] * 3 - 3] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[pole[i][1] * 3 - 2][pole[i][1] * 3 - 2] = (cos_y ** 2 * model[i] * area[i] / longs[
                i] + cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[pole[i][1] * 3 - 3][pole[i][1] * 3 - 1] = (6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[pole[i][1] * 3 - 1][pole[i][1] * 3 - 3] = (6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[pole[i][1] * 3 - 1][pole[i][1] * 3 - 1] = (4 * model[i] * inertia[i] / longs[i])
            k[pole[i][1] * 3 - 2][pole[i][1] * 3 - 1] = (-6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[pole[i][1] * 3 - 1][pole[i][1] * 3 - 2] = (-6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            stiff.append(k)
        g_stiff = Structures.zeros(3*n_joints, 3*n_joints)
        for bars in range(len(pole)):
            g_stiff = Structures.sum_matrix(g_stiff, stiff[bars])
        return g_stiff

    @staticmethod
    def local_stiffness(joints, pole):
        l_stiff = []
        longs = Structures.longitude(joints, pole)
        model = [i[4] for i in pole]
        inertia = Structures.moment_inertia(pole)
        area = Structures.area_poles(pole)
        for i in range(len(pole)):
            cos_x = (joints[pole[i][1] - 1][0] - joints[pole[i][0] - 1][0]) / longs[i]
            cos_y = (joints[pole[i][1] - 1][1] - joints[pole[i][0] - 1][1]) / longs[i]
            k = Structures.zeros(6, 6)
            k[0][0] = (cos_x ** 2 * model[i] * area[i] / longs[i] +
                       cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[0][1] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[i] /
                        longs[i] ** 3) * cos_y * cos_x)
            k[1][0] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[i] /
                        longs[i] ** 3) * cos_y * cos_x)
            k[1][1] = (cos_y ** 2 * model[i] * area[i] / longs[i] +
                       cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[0][2] = (-6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[2][0] = (-6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[2][2] = (4 * model[i] * inertia[i] / longs[i])
            k[1][2] = (6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[2][1] = (6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[0][3] = - (cos_x ** 2 * model[i] * area[i] / longs[
                i] + cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[0][4] = - ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[1][3] = - ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[1][4] = - (cos_y ** 2 * model[i] * area[i] / longs[
                i] + cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[0][5] = (-6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[2][3] = (6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[1][5] = (6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[2][4] = (-6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[2][5] = (2 * model[i] * inertia[i] / longs[i])
            k[3][0] = - (cos_x ** 2 * model[i] * area[i] / longs[
                i] + cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[4][0] = - ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[3][1] = - ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[4][1] = - (cos_y ** 2 * model[i] * area[i] / longs[
                i] + cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[5][0] = (-6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[3][2] = (6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[5][1] = (6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[4][2] = (-6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[5][2] = (2 * model[i] * inertia[i] / longs[i])
            k[3][3] = (cos_x ** 2 * model[i] * area[i] / longs[
                i] + cos_y ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[3][4] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[4][3] = ((model[i] * area[i] / longs[i] - 12 * model[i] * inertia[
                i] / longs[i] ** 3) * cos_y * cos_x)
            k[4][4] = (cos_y ** 2 * model[i] * area[i] / longs[
                i] + cos_x ** 2 * 12 * model[i] * inertia[i] / longs[i] ** 3)
            k[3][5] = (6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[5][3] = (6 * cos_y * model[i] * inertia[i] / longs[i] ** 2)
            k[5][5] = (4 * model[i] * inertia[i] / longs[i])
            k[4][5] = (-6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            k[5][4] = (-6 * cos_x * model[i] * inertia[i] / longs[i] ** 2)
            l_stiff.append(k)
        return l_stiff

    @staticmethod
    def forces_by_internal_loads(joints, pole, internal_loads):
        forces_internal_load = []
        longs = Structures.longitude(joints, pole)
        n_joints = len(joints)
        for i in range(len(internal_loads)):
            n_ele = internal_loads[i][0]-1
            forces_by_punctual = Structures.zeros(3*n_joints, 1)
            forces_by_distributed = Structures.zeros(3*n_joints, 1)
            forces_by_momentum = Structures.zeros(3*n_joints, 1)
            condition_joint = internal_loads[i][1]
            long = longs[n_ele]
            cos_x = (joints[pole[i][1] - 1][0] - joints[pole[i][0] - 1][0]) / longs[i]
            sen_x = (joints[pole[i][1] - 1][1] - joints[pole[i][0] - 1][1]) / longs[i]
            a = 0
            b = 1
            c = 1
            if (cos_x < 0) and (sen_x < 0 or sen_x > 0):
                a = 1
                b = 0
                c = -1
            elif (sen_x < 0 and cos_x == 0) or (sen_x == 0 and cos_x < 0):
                a = 1
                b = 0
                c = -1
            if condition_joint == 'punctual':
                force = internal_loads[i][3]
                x_p = internal_loads[i][2]
                forces_by_punctual[pole[n_ele][a] * 3 - 3][
                    0] = -c*sen_x * force * (long - x_p) ** 2 * (3 - 2 * (long - x_p) / long) / long ** 2
                forces_by_punctual[pole[n_ele][a] * 3 - 2][
                    0] = c*cos_x * force * (long - x_p) ** 2 * (3 - 2 * (long - x_p) / long) / long ** 2
                forces_by_punctual[pole[n_ele][a] * 3 - 1][
                    0] = force * x_p * (long - x_p) ** 2 / long ** 2
                forces_by_punctual[pole[n_ele][b] * 3 - 3][
                    0] = - c*sen_x * force * x_p ** 2 * (3 - 2 * x_p / long) / long ** 2
                forces_by_punctual[pole[n_ele][b] * 3 - 2][
                    0] = c*cos_x * force * x_p ** 2 * (3 - 2 * x_p / long) / long ** 2
                forces_by_punctual[pole[n_ele][b] * 3 - 1][
                    0] = -force * x_p ** 2 * (long - x_p) / long ** 2

                forces_internal_load.append(forces_by_punctual)
            if condition_joint == 'distributed':
                force = internal_loads[i][4]
                x_a = internal_loads[i][2]
                x_b = internal_loads[i][3]
                x_c = x_a + x_b/2
                forces_by_distributed[pole[n_ele][a] * 3 - 3][
                    0] = -c*sen_x*force*x_b*(1-3*(x_c/long)**2-0.25*(x_b**2)/(long**2)+2*x_c*(x_c**2+0.25*x_b**2)/long**3)
                forces_by_distributed[pole[n_ele][a] * 3 - 2][
                    0] = c*cos_x*force*x_b*(1-3*(x_c/long)**2-0.25*(x_b**2)/(long**2)+2*x_c*(x_c**2+0.25*x_b**2)/long**3)
                forces_by_distributed[pole[n_ele][a] * 3 - 1][
                    0] = force * x_b * (x_c*(long-x_c)**2+x_b**2*(long-3*(long-x_c))/12) / long ** 2
                forces_by_distributed[pole[n_ele][b] * 3 - 3][
                    0] = -c*sen_x*force*x_b*(3*(x_c/long)**2+0.25*(x_b**2)/(long**2)-2*x_c*(x_c**2+0.25*x_b**2)/long**3)
                forces_by_distributed[pole[n_ele][b] * 3 - 2][
                    0] = c*cos_x*force*x_b*(3*(x_c/long)**2+0.25*(x_b**2)/(long**2)-2*x_c*(x_c**2+0.25*x_b**2)/long**3)
                forces_by_distributed[pole[n_ele][b] * 3 - 1][
                    0] = -force * x_b * (x_c**2*(long-x_c)+x_b**2*(long-3*x_c)/12) / long ** 2
                forces_internal_load.append(forces_by_distributed)
            if condition_joint == 'momentum':
                momentum = internal_loads[i][3]
                x_m = internal_loads[i][2]
                forces_by_momentum[pole[n_ele][a] * 3 - 3][
                    0] = -c*sen_x * 6 * momentum * x_m * (long - x_m) / long ** 3
                forces_by_momentum[pole[n_ele][a] * 3 - 2][
                    0] = c*cos_x * 6 * momentum * x_m * (long - x_m) / long ** 3
                forces_by_momentum[pole[n_ele][a] * 3 - 1][
                    0] = momentum * (long - x_m) * (2 - 3 * ((long-x_m)/long)) / long
                forces_by_momentum[pole[n_ele][b] * 3 - 3][
                    0] = c*sen_x * 6 * momentum * x_m * (long - x_m) / long ** 3
                forces_by_momentum[pole[n_ele][b] * 3 - 2][
                    0] = - c*cos_x * 6 * momentum * x_m * (long - x_m) / long ** 3
                forces_by_momentum[pole[n_ele][b] * 3 - 1][
                    0] = momentum * x_m * (2 - 3*x_m/long) / long
                forces_internal_load.append(forces_by_momentum)
        global_forces_internal_load = Structures.zeros(3*n_joints, 1)
        for i in range(len(forces_internal_load)):
            global_forces_internal_load = Structures.sum_matrix(global_forces_internal_load, forces_internal_load[i])
        return global_forces_internal_load

    @staticmethod
    def forces_by_external_loads(external_loads):
        forces_external_load = []
        n_joints = len(external_loads)
        for i in range(n_joints):
            forces_external_load.append([external_loads[i][0]])
            forces_external_load.append([external_loads[i][1]])
            forces_external_load.append([external_loads[i][2]])
        return forces_external_load

    @staticmethod
    def grades_for_freedom(joints):
        """
        |:param joints: 
        |    List of coordinates and condition of the board:
        |    0: simple
        |    1: mobile
        |    2: fixed
        |    3: recessed 
        |:return: g_f
        |    List of grades of freedom of the structure
        """
        g_f = []
        count = 1
        for i in joints:
            if i[2] == 1:
                g_f.append(count*3-1)
            elif i[2] == 2:
                g_f.append(count*3-2)
                g_f.append(count*3-1)
            elif i[2] == 3:
                g_f.append(count*3-2)
                g_f.append(count*3-1)
                g_f.append(count*3)
            count += 1
        return g_f

    @staticmethod
    def displacement_nods(joints, pole, internal_loads, external_loads):
        grades_freedom = Structures.grades_for_freedom(joints)
        ngl = len(grades_freedom)
        m_stiffness = Structures.global_stiffness(joints, pole)
        f_internal = Structures.forces_by_internal_loads(joints, pole, internal_loads)
        f_external = Structures.forces_by_external_loads(external_loads)
        f_t = Structures.sum_matrix(f_external, Structures.scalar_product(f_internal, -1))
        for i in range(ngl):
            Structures.remove_column(m_stiffness, grades_freedom[ngl - i - 1])
            Structures.remove_row(m_stiffness, grades_freedom[ngl - i - 1])
            Structures.remove_row(f_t, grades_freedom[ngl - i - 1])
        displacements = Structures.matrix_solve(m_stiffness, f_t)
        for i in grades_freedom:
            displacements.insert(i-1, [0])
        return displacements


joint = [[0.0, 4.0, 0], [3.0, 4.0, 0], [0.0, 0.0, 3], [3.0, 0.0, 3]]
poles = [[1, 2, 1, 1, 1], [3, 1, 1, 1, 1], [2, 4, 1, 1, 1]]
internal_load = [[1, 'distributed', 0, 3, 1.2], [2, 'punctual', 2, 2], [3, 'momentum', 1.5, 1]]
external_load = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
"""print(Structures.forces_by_internal_loads(joint, poles, internal_load))
g_s = Structures.global_stiffness(joint, poles)
print(g_s)
f_e = Structures.forces_by_internal_loads(joint, poles, internal_load)

help(Structures)"""
print(Structures.displacement_nods(joint, poles, internal_load, external_load))
