import numpy as np
cimport cython
from numpy cimport ndarray
cimport numpy as np
from libc.math cimport sqrt, exp
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cdef double part_distance(double dx, double dy, double dz):
    cdef double r = sqrt(dx * dx + dy * dy + dz * dz)
    return r


cdef double smooth_distance(double dx, double dy, double dz,
                            double eps_smooth):
    cdef double delta_soft = sqrt(dx * dx + dy * dy + dz * dz
                                  + eps_smooth * eps_smooth)
    return delta_soft


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef int_Ps_to_P(np.ndarray[DTYPE_t, ndim=2] Particles,
                int n1, int n2,
                int n_1, int n_2, double eps_smooth):
    assert Particles.dtype == DTYPE
    cdef double a_x
    cdef double a_y
    cdef double a_z
    cdef double phi
    cdef double r_1
    cdef double r_3
    cdef double m
    cdef double dx0
    cdef double x1
    cdef double dy0
    cdef double y1
    cdef double dz0
    cdef double z1
    cdef int part_num
    cdef int num
    cdef int counter = 0
    cdef ndarray[DTYPE_t, ndim=2] A = np.zeros([n2 - n1, 4], dtype=DTYPE)
    for part_num in range(n1, n2):
        a_x = 0
        a_y = 0
        a_z = 0
        phi = 0
        if (part_num >= n_1) and (part_num < n_2):
            for num in range(n_1, n_2):
                if not num == part_num:
                    m = Particles[num, 6]
                    dx0 = Particles[part_num, 0]
                    x1 = Particles[num, 0]
                    dy0 = Particles[part_num, 1]
                    y1 = Particles[num, 1]
                    dz0 = Particles[part_num, 2]
                    z1 = Particles[num, 2]
                    dx0 = x1 - dx0
                    dy0 = y1 - dy0
                    dz0 = z1 - dz0
                    r_1 = smooth_distance(dx0, dy0, dz0, eps_smooth)
                    r_3 = m / (r_1 * r_1 * r_1)
                    a_x += dx0 * r_3
                    a_y += dy0 * r_3
                    a_z += dz0 * r_3
                    phi += m / r_1
        else:
            for num in range(n_1, n_2):
                m = Particles[num, 6]
                dx0 = Particles[part_num, 0]
                x1 = Particles[num, 0]
                dy0 = Particles[part_num, 1]
                y1 = Particles[num, 1]
                dz0 = Particles[part_num, 2]
                z1 = Particles[num, 2]
                dx0 = x1 - dx0
                dy0 = y1 - dy0
                dz0 = z1 - dz0
                r_1 = smooth_distance(dx0, dy0, dz0, eps_smooth)
                r_3 = m / (r_1 * r_1 * r_1)
                a_x += dx0 * r_3
                a_y += dy0 * r_3
                a_z += dz0 * r_3
                phi += m / r_1
        A[counter, 0] = a_x
        A[counter, 1] = a_y
        A[counter, 2] = a_z
        A[counter, 3] = - phi
        counter += 1
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef int_C_to_P(np.ndarray[DTYPE_t, ndim=2] Particles,
               np.ndarray[DTYPE_t, ndim=2] Mass_center,
               int Part_num, int cell_num):
    assert Particles.dtype == DTYPE and Mass_center.dtype == DTYPE
    # Функция, рассчитывающая ускорение частицы под номером Part_num,
    # полученное за счет гравитационного мультипольного взаимодействия с
    # частицами в ячейке с номером cell_num.
    cdef ndarray[DTYPE_t, ndim=1] A = np.zeros([4], dtype=DTYPE)
    cdef double m = Mass_center[cell_num, 6]
    cdef double dx0 = Particles[Part_num, 0]
    cdef double x1 = Mass_center[cell_num, 0]
    cdef double dy0 = Particles[Part_num, 1]
    cdef double y1 = Mass_center[cell_num, 1]
    cdef double dz0 = Particles[Part_num, 2]
    cdef double z1 = Mass_center[cell_num, 2]
    dx0 = x1 - dx0
    dy0 = y1 - dy0
    dz0 = z1 - dz0
    cdef double r_1 = part_distance(dx0, dy0, dz0)
    cdef double r_3 = m / (r_1 * r_1 * r_1)
    dx0 = dx0 * r_3
    dy0 = dy0 * r_3
    dz0 = dz0 * r_3
    cdef double phi = - m / r_1
#    cell_to_body += quadrupole(Mass_center, cell_num, r_1, r_3,
#                               delta_x, delta_y, delta_z)
    A[0] = dx0
    A[1] = dy0
    A[2] = dz0
    A[3] = phi
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef g_force_Newton(np.ndarray[DTYPE_t, ndim=2] Particles,
                     int part_num, int total_part, double smooth):
    assert Particles.dtype == DTYPE
    # Ускорение по Ньютону
    cdef ndarray[DTYPE_t, ndim=2] A = np.zeros([total_part, 4], dtype=DTYPE)
    cdef double m
    cdef double m1
    cdef double dx0
    cdef double x1
    cdef double dy0
    cdef double y1
    cdef double dz0
    cdef double z1
    cdef double r_1
    cdef double r_3
    cdef double r3
    cdef int k
    cdef int l
    for k in range(part_num):
        for l in range(total_part):
            if not l == k:
                m = Particles[k, 6]
                dx0 = Particles[l, 0]
                x1 = Particles[k, 0]
                dy0 = Particles[l, 1]
                y1 = Particles[k, 1]
                dz0 = Particles[l, 2]
                z1 = Particles[k, 2]
                dx0 = x1 - dx0
                dy0 = y1 - dy0
                dz0 = z1 - dz0
                r_1 = smooth_distance(dx0, dy0, dz0, smooth)
                r_3 = 1 / (r_1 * r_1 * r_1)
                if l >= part_num:
                    m1 = Particles[l, 6]
                    r3 = r_3 * m1
                    A[k, 0] += - dx0 * r3
                    A[k, 1] += - dy0 * r3
                    A[k, 2] += - dz0 * r3
                    A[k, 3] += - m1 / r_1
                r_3 = r_3 * m
                A[l, 0] += dx0 * r_3
                A[l, 1] += dy0 * r_3
                A[l, 2] += dz0 * r_3
                A[l, 3] += - m / r_1
    return A

@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef interaction_Ps2P_plus_Yukawa(np.ndarray[DTYPE_t, ndim=2] Particles,
                                   int n1, int n2,
                                   int n_1, int n_2, double eps_smooth,
                                   double alpha, double lambd):
    assert Particles.dtype == DTYPE
    cdef double a_x
    cdef double a_y
    cdef double a_z
    cdef double phi
    cdef double r_1
    cdef double r1
    cdef double r_2
    cdef double m
    cdef double dx0
    cdef double x1
    cdef double dy0
    cdef double y1
    cdef double dz0
    cdef double z1
    cdef double exp_add
    cdef int part_num
    cdef int num
    cdef int counter = 0
    cdef ndarray[DTYPE_t, ndim=2] A = np.zeros([n2 - n1, 4], dtype=DTYPE)
    for part_num in range(n1, n2):
        a_x = 0
        a_y = 0
        a_z = 0
        phi = 0
        if (part_num >= n_1) and (part_num < n_2):
            for num in range(n_1, n_2):
                if not num == part_num:
                    m = Particles[num, 6]
                    dx0 = Particles[part_num, 0]
                    x1 = Particles[num, 0]
                    dy0 = Particles[part_num, 1]
                    y1 = Particles[num, 1]
                    dz0 = Particles[part_num, 2]
                    z1 = Particles[num, 2]
                    dx0 = x1 - dx0
                    dy0 = y1 - dy0
                    dz0 = z1 - dz0
                    r1 = smooth_distance(dx0, dy0, dz0, eps_smooth)
                    r_1 = 1 / r1
                    r_2 = m * r_1 * r_1
                    exp_add = alpha * exp(- r1 / lambd)
                    r_2 = r_2 * (r_1 + exp_add * (r_1 + (1 / lambd)))
                    a_x += dx0 * r_2
                    a_y += dy0 * r_2
                    a_z += dz0 * r_2
                    phi += (m * r_1) * (1 + exp_add)
        else:
            for num in range(n_1, n_2):
                m = Particles[num, 6]
                dx0 = Particles[part_num, 0]
                x1 = Particles[num, 0]
                dy0 = Particles[part_num, 1]
                y1 = Particles[num, 1]
                dz0 = Particles[part_num, 2]
                z1 = Particles[num, 2]
                dx0 = x1 - dx0
                dy0 = y1 - dy0
                dz0 = z1 - dz0
                r1 = smooth_distance(dx0, dy0, dz0, eps_smooth)
                r_1 = 1 / r1
                r_2 = m * r_1 * r_1
                exp_add = alpha * exp(- r1 / lambd)
                r_2 = r_2 * (r_1 + exp_add * (r_1 + (1 / lambd)))
                a_x += dx0 * r_2
                a_y += dy0 * r_2
                a_z += dz0 * r_2
                phi += (m * r_1) * (1 + exp_add)
        A[counter, 0] = a_x
        A[counter, 1] = a_y
        A[counter, 2] = a_z
        A[counter, 3] = - phi
        counter += 1
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef interaction_C2P_plus_Yukawa(np.ndarray[DTYPE_t, ndim=2] Particles,
                                  np.ndarray[DTYPE_t, ndim=2] Mass_center,
                                  double alpha, double lambd,
                                  int Part_num, int cell_num):
    assert Particles.dtype == DTYPE and Mass_center.dtype == DTYPE
    cdef ndarray[DTYPE_t, ndim=1] A = np.zeros([4], dtype=DTYPE)
    cdef double m = Mass_center[cell_num, 6]
    cdef double dx0 = Particles[Part_num, 0]
    cdef double x1 = Mass_center[cell_num, 0]
    cdef double dy0 = Particles[Part_num, 1]
    cdef double y1 = Mass_center[cell_num, 1]
    cdef double dz0 = Particles[Part_num, 2]
    cdef double z1 = Mass_center[cell_num, 2]
    cdef double exp_add
    cdef double r_2
    dx0 = x1 - dx0
    dy0 = y1 - dy0
    dz0 = z1 - dz0
    cdef double r1 = part_distance(dx0, dy0, dz0)
    cdef double r_1 = 1 / r1
    exp_add = alpha * exp(- r1 / lambd)
    r_2 = m * r_1 * r_1 * (r_1 + exp_add * (r_1 + (1 / lambd)))
    dx0 = dx0 * r_2
    dy0 = dy0 * r_2
    dz0 = dz0 * r_2
    cdef double phi = (m * r_1) * (1 + exp_add)
    A[0] = dx0
    A[1] = dy0
    A[2] = dz0
    A[3] = phi
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef interaction_Ps2P_Brance_Dicke(np.ndarray[DTYPE_t, ndim=2] Particles,
                                    int n1, int n2,
                                    int n_1, int n_2, double eps_smooth,
                                    double gamma):
    assert Particles.dtype == DTYPE
    cdef double a_x
    cdef double a_y
    cdef double a_z
    cdef double phi
    cdef double r_1
    cdef double r1
    cdef double r_3
    cdef double m
    cdef double dx0
    cdef double x1
    cdef double dy0
    cdef double y1
    cdef double dz0
    cdef double z1
    cdef int part_num
    cdef int num
    cdef int counter = 0
    cdef ndarray[DTYPE_t, ndim=2] A = np.zeros([n2 - n1, 4], dtype=DTYPE)
    for part_num in range(n1, n2):
        a_x = 0
        a_y = 0
        a_z = 0
        phi = 0
        if (part_num >= n_1) and (part_num < n_2):
            for num in range(n_1, n_2):
                if not num == part_num:
                    m = Particles[num, 6]
                    dx0 = Particles[part_num, 0]
                    x1 = Particles[num, 0]
                    dy0 = Particles[part_num, 1]
                    y1 = Particles[num, 1]
                    dz0 = Particles[part_num, 2]
                    z1 = Particles[num, 2]
                    dx0 = x1 - dx0
                    dy0 = y1 - dy0
                    dz0 = z1 - dz0
                    r1 = smooth_distance(dx0, dy0, dz0, eps_smooth)
                    r_1 = 1 / r1
                    r_3 = m * r_1 * r_1 * r_1
                    a_x += dx0 * r_3 * (1 + 2 * gamma * m * r_1)
                    a_y += dy0 * r_3 * (1 + 2 * gamma * m * r_1)
                    a_z += dz0 * r_3 * (1 + 2 * gamma * m * r_1)
                    phi += (m * r_1) * (1 + gamma * m * r_1)
        else:
            for num in range(n_1, n_2):
                m = Particles[num, 6]
                dx0 = Particles[part_num, 0]
                x1 = Particles[num, 0]
                dy0 = Particles[part_num, 1]
                y1 = Particles[num, 1]
                dz0 = Particles[part_num, 2]
                z1 = Particles[num, 2]
                dx0 = x1 - dx0
                dy0 = y1 - dy0
                dz0 = z1 - dz0
                r1 = smooth_distance(dx0, dy0, dz0, eps_smooth)
                r_1 = 1 / r1
                r_3 = m * r_1 * r_1 * r_1
                a_x += dx0 * r_3 * (1 + 2 * gamma * m * r_1)
                a_y += dy0 * r_3 * (1 + 2 * gamma * m * r_1)
                a_z += dz0 * r_3 * (1 + 2 * gamma * m * r_1)
                phi += (m * r_1) * (1 + gamma * m * r_1)
        A[counter, 0] = a_x
        A[counter, 1] = a_y
        A[counter, 2] = a_z
        A[counter, 3] = - phi
        counter += 1
    return A

@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef interaction_C2P_Brance_Dicke(np.ndarray[DTYPE_t, ndim=2] Particles,
                                   np.ndarray[DTYPE_t, ndim=2] Mass_center,
                                   double gamma,
                                   int Part_num, int cell_num):
    assert Particles.dtype == DTYPE and Mass_center.dtype == DTYPE
    cdef ndarray[DTYPE_t, ndim=1] A = np.zeros([4], dtype=DTYPE)
    cdef double m = Mass_center[cell_num, 6]
    cdef double m_sum = Mass_center[cell_num, 23]
    cdef double dx0 = Particles[Part_num, 0]
    cdef double x1 = Mass_center[cell_num, 0]
    cdef double dy0 = Particles[Part_num, 1]
    cdef double y1 = Mass_center[cell_num, 1]
    cdef double dz0 = Particles[Part_num, 2]
    cdef double z1 = Mass_center[cell_num, 2]
    cdef double r_3
    dx0 = x1 - dx0
    dy0 = y1 - dy0
    dz0 = z1 - dz0
    cdef double r1 = part_distance(dx0, dy0, dz0)
    cdef double r_1 = 1 / r1
    r_3 = r_1 * r_1 * r_1
    dx0 = dx0 * r_3 * (m + 2 * gamma * m_sum * r_1)
    dy0 = dy0 * r_3 * (m + 2 * gamma * m_sum * r_1)
    dz0 = dz0 * r_3 * (m + 2 * gamma * m_sum * r_1)
    cdef double phi = r_1 * (m + gamma * m_sum * r_1)
    A[0] = dx0
    A[1] = dy0
    A[2] = dz0
    A[3] = phi
    return A


# def quadrupole(Mass_center, num, r_1, r_3, delta_x, delta_y, delta_z):
#    # Функция, расчитывающая квадрупольный вклад
#    r_5 = r_3 * r_1 * r_1
#    r_7 = r_5 * r_1 * r_1
#    DR = (Mass_center[num, 7] * delta_x * delta_y
#          + Mass_center[num, 8] * delta_x * delta_z
#          + Mass_center[num, 9] * delta_y * delta_z) * 5
#    a_x = - (Mass_center[num, 7] * delta_y + Mass_center[num, 8] * delta_z) \
#        / r_5 + DR * delta_x / r_7
#    a_y = - (Mass_center[num, 7] * delta_x + Mass_center[num, 9] * delta_z) \
#        / r_5 + DR * delta_y / r_7
#    a_z = - (Mass_center[num, 8] * delta_x + Mass_center[num, 9] * delta_y) \
#        / r_5 + DR * delta_z / r_7
#    phi = DR / (5 * r_5)
#    return np.array([a_x, a_y, a_z, - phi])
