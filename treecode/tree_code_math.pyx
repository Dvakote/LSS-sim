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

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef sqr_dist(np.ndarray[DTYPE_t, ndim=2] Mass_center,
               int current_cell, int cell_num):
#    assert Mass_center.dtype == DTYPE
    cdef:
        double r_2
        double x1 = Mass_center[current_cell, 3]
        double x0 = Mass_center[cell_num, 3]
        double y1 = Mass_center[current_cell, 4]
        double y0 = Mass_center[cell_num, 4]
        double z1 = Mass_center[current_cell, 5]
        double z0 = Mass_center[cell_num, 5]
    x1 = x1 - x0
    y1 = y1 - y0
    z1 = z1 - z0
    r_2 = x1 * x1 + y1 * y1 + z1 * z1
    return r_2


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef Ps2Cell_fast(np.ndarray[DTYPE_t, ndim=2] Y,
                   double CELL_LENGTH,
                   int Y_size, int order_n,
                   int n, int n_max):
#    assert Y.dtype == DTYPE
    cdef np.ndarray[DTYPE_t, ndim=2] R_local = np.zeros([n_max, 18], dtype=DTYPE)
    cdef:
        int cell_num
        int cell_x
        int n1
        int n2
        int part_num = 0
        int part_count = 0
        int Particle_in_cell_num = int(Y[part_num, 11]+0.001)
        double x
        double y
        double z
        double m
        double x_sum
        double y_sum
        double z_sum
        double m_sum
        double m_sqr_sum
        double x_center = 0
        double y_center = 0
        double z_center = 0
        double L_2 = (9 - 0.000001) * CELL_LENGTH * CELL_LENGTH
    while Particle_in_cell_num < 0:
        part_count += 1
        part_num += 1
        if part_num == Y_size:
            break
        Particle_in_cell_num = int(Y[part_num, 11]+0.001)
    for cell_num in range(n_max):
        x_sum = 0
        y_sum = 0
        z_sum = 0
        m_sum = 0
        m_sqr_sum = 0
        if not part_num == Y_size:
            while Particle_in_cell_num == cell_num:
                x = Y[part_num, 0]
                y = Y[part_num, 1]
                z = Y[part_num, 2]
                m = Y[part_num, 6]
                x_sum += x * m
                y_sum += y * m
                z_sum += z * m
                m_sum += m
                m_sqr_sum += m * m
                part_num += 1
                if part_num == Y_size:
                    break
                Particle_in_cell_num = int(Y[part_num, 11]+0.001)
        n1 = part_count
        n2 = part_num
        part_count = part_num
        if not m_sum == 0:
            # Расчет положения центра масс ячейки
            x_sum = x_sum / m_sum
            y_sum = y_sum / m_sum
            z_sum = z_sum / m_sum
            # Расчет положения геометрического центра ячейки
            cell_x = cell_num // (n * n)
            x_center = CELL_LENGTH * (0.5 + cell_x)
            y_center = CELL_LENGTH * (0.5 + ((cell_num // n) - cell_x * n))
            z_center = CELL_LENGTH * (0.5 + (cell_num % n))
        R_local[cell_num, 0] = x_sum
        R_local[cell_num, 1] = y_sum
        R_local[cell_num, 2] = z_sum
        R_local[cell_num, 3] = x_center
        R_local[cell_num, 4] = y_center
        R_local[cell_num, 5] = z_center
        R_local[cell_num, 6] = m_sum
        R_local[cell_num, 7] = m_sqr_sum
        R_local[cell_num, 8] = L_2
        R_local[cell_num, 9] = order_n
        R_local[cell_num, 10] = n1
        R_local[cell_num, 11] = n2
    return R_local


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef int_Ps_to_P(np.ndarray[DTYPE_t, ndim=2] Particles,
                int n1, int n2,
                int n_1, int n_2, double eps_smooth):
#    assert Particles.dtype == DTYPE
    cdef:
        double a_x
        double a_y
        double a_z
        double phi
        double r_1
        double r_3
        double m
        double dx0
        double x1
        double dy0
        double y1
        double dz0
        double z1
        int part_num
        int num
        int counter = 0
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.zeros([n2 - n1, 4], dtype=DTYPE)
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
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef int_C_to_P(np.ndarray[DTYPE_t, ndim=2, mode='c'] Particles,
               np.ndarray[DTYPE_t, ndim=2, mode='c'] Mass_center,
               int Part_num, int cell_num):
#    assert Particles.dtype == DTYPE and Mass_center.dtype == DTYPE
    # Функция, рассчитывающая ускорение частицы под номером Part_num,
    # полученное за счет гравитационного мультипольного взаимодействия с
    # частицами в ячейке с номером cell_num.
    cdef np.ndarray[DTYPE_t, ndim=1] A = np.zeros([4], dtype=DTYPE)
    cdef:
        double m = Mass_center[cell_num, 6]
        double dx0 = Particles[Part_num, 0]
        double x1 = Mass_center[cell_num, 0]
        double dy0 = Particles[Part_num, 1]
        double y1 = Mass_center[cell_num, 1]
        double dz0 = Particles[Part_num, 2]
        double z1 = Mass_center[cell_num, 2]
        double r_1
        double r_3
        double phi
    dx0 = x1 - dx0
    dy0 = y1 - dy0
    dz0 = z1 - dz0
    r_1 = part_distance(dx0, dy0, dz0)
    r_3 = m / (r_1 * r_1 * r_1)
    dx0 = dx0 * r_3
    dy0 = dy0 * r_3
    dz0 = dz0 * r_3
    phi = - m / r_1
#    cell_to_body += quadrupole(Mass_center, cell_num, r_1, r_3,
#                               delta_x, delta_y, delta_z)
    A[0] = dx0
    A[1] = dy0
    A[2] = dz0
    A[3] = phi
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef g_force_Newton(np.ndarray[DTYPE_t, ndim=2] Particles,
                     int part_num, int total_part, double smooth):
#    assert Particles.dtype == DTYPE
    # Ускорение по Ньютону
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.zeros([total_part, 4], dtype=DTYPE)
    cdef:
        double m
        double m1
        double dx0
        double x1
        double dy0
        double y1
        double dz0
        double z1
        double r_1
        double r_3
        double r3
        int k
        int l
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
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef interaction_Ps2P_plus_Yukawa(np.ndarray[DTYPE_t, ndim=2] Particles,
                                   int n1, int n2,
                                   int n_1, int n_2, double eps_smooth,
                                   double alpha, double lambd):
#    assert Particles.dtype == DTYPE
    cdef:
        double a_x
        double a_y
        double a_z
        double phi
        double r_1
        double r1
        double r_2
        double m
        double dx0
        double x1
        double dy0
        double y1
        double dz0
        double z1
        double exp_add
        int part_num
        int num
        int counter = 0
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.zeros([n2 - n1, 4], dtype=DTYPE)
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
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef interaction_C2P_plus_Yukawa(np.ndarray[DTYPE_t, ndim=2] Particles,
                                  np.ndarray[DTYPE_t, ndim=2] Mass_center,
                                  double alpha, double lambd,
                                  int Part_num, int cell_num):
#    assert Particles.dtype == DTYPE and Mass_center.dtype == DTYPE
    cdef np.ndarray[DTYPE_t, ndim=1, mode='c'] A = np.zeros([4], dtype=DTYPE)
    cdef:
        double m = Mass_center[cell_num, 6]
        double dx0 = Particles[Part_num, 0]
        double x1 = Mass_center[cell_num, 0]
        double dy0 = Particles[Part_num, 1]
        double y1 = Mass_center[cell_num, 1]
        double dz0 = Particles[Part_num, 2]
        double z1 = Mass_center[cell_num, 2]
        double exp_add
        double r_2
        double r1
        double r_1
        double phi
    dx0 = x1 - dx0
    dy0 = y1 - dy0
    dz0 = z1 - dz0
    r1 = part_distance(dx0, dy0, dz0)
    r_1 = 1 / r1
    exp_add = alpha * exp(- r1 / lambd)
    r_2 = m * r_1 * r_1 * (r_1 + exp_add * (r_1 + (1 / lambd)))
    dx0 = dx0 * r_2
    dy0 = dy0 * r_2
    dz0 = dz0 * r_2
    phi = (m * r_1) * (1 + exp_add)
    A[0] = dx0
    A[1] = dy0
    A[2] = dz0
    A[3] = phi
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef g_force_plusYukawa(np.ndarray[DTYPE_t, ndim=2] Particles,
                         int part_num, int total_part,
                         double smooth, double alpha, double lambd):
#    assert Particles.dtype == DTYPE
    # Ускорение по Ньютону
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.zeros([total_part, 4], dtype=DTYPE)
    cdef:
        double m
        double m1
        double dx0
        double x1
        double dy0
        double y1
        double dz0
        double z1
        double r1
        double r_1
        double r_2
        double r2
        double exp_add
        int k
        int l
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
                r1 = smooth_distance(dx0, dy0, dz0, smooth)
                r_1 = 1 / r1
                r_2 = r_1 * r_1
                exp_add = alpha * exp(- r1 / lambd)
                r_2 = r_2 * (r_1 + exp_add * (r_1 + (1 / lambd)))
                if l >= part_num:
                    m1 = Particles[l, 6]
                    r2 = r_2 * m1
                    A[k, 0] += - dx0 * r2
                    A[k, 1] += - dy0 * r2
                    A[k, 2] += - dz0 * r2
                    A[k, 3] += - (m1 * r_1) * (1 + exp_add)
                r_2 = r_2 * m
                A[l, 0] += dx0 * r_2
                A[l, 1] += dy0 * r_2
                A[l, 2] += dz0 * r_2
                A[l, 3] += - (m * r_1) * (1 + exp_add)
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef interaction_Ps2P_Brance_Dicke(np.ndarray[DTYPE_t, ndim=2] Particles,
                                    int n1, int n2,
                                    int n_1, int n_2, double eps_smooth,
                                    double gamma):
#    assert Particles.dtype == DTYPE
    cdef:
        double a_x
        double a_y
        double a_z
        double phi
        double r_1
        double r1
        double r_3
        double m
        double dx0
        double x1
        double dy0
        double y1
        double dz0
        double z1
        int part_num
        int num
        int counter = 0
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.zeros([n2 - n1, 4], dtype=DTYPE)
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
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef interaction_C2P_Brance_Dicke(np.ndarray[DTYPE_t, ndim=2] Particles,
                                   np.ndarray[DTYPE_t, ndim=2] Mass_center,
                                   double gamma,
                                   int Part_num, int cell_num):
#    assert Particles.dtype == DTYPE and Mass_center.dtype == DTYPE
    cdef np.ndarray[DTYPE_t, ndim=1, mode='c'] A = np.zeros([4], dtype=DTYPE)
    cdef:
        double m = Mass_center[cell_num, 6]
        double m_sum = Mass_center[cell_num, 7]
        double dx0 = Particles[Part_num, 0]
        double x1 = Mass_center[cell_num, 0]
        double dy0 = Particles[Part_num, 1]
        double y1 = Mass_center[cell_num, 1]
        double dz0 = Particles[Part_num, 2]
        double z1 = Mass_center[cell_num, 2]
        double r1
        double r_1
        double r_3
        double phi
    dx0 = x1 - dx0
    dy0 = y1 - dy0
    dz0 = z1 - dz0
    r1 = part_distance(dx0, dy0, dz0)
    r_1 = 1 / r1
    r_3 = r_1 * r_1 * r_1
    dx0 = dx0 * r_3 * (m + 2 * gamma * m_sum * r_1)
    dy0 = dy0 * r_3 * (m + 2 * gamma * m_sum * r_1)
    dz0 = dz0 * r_3 * (m + 2 * gamma * m_sum * r_1)
    phi = r_1 * (m + gamma * m_sum * r_1)
    A[0] = dx0
    A[1] = dy0
    A[2] = dz0
    A[3] = phi
    return A


@cython.cdivision(True)
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef g_force_BD(np.ndarray[DTYPE_t, ndim=2] Particles,
                 int part_num, int total_part,
                 double smooth, double gamma):
#    assert Particles.dtype == DTYPE
    # Ускорение по Ньютону
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.zeros([total_part, 4], dtype=DTYPE)
    cdef:
        double m
        double m1
        double dx0
        double x1
        double dy0
        double y1
        double dz0
        double z1
        double r1
        double r_1
        double r_3
        double r3
        int k
        int l
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
                r1 = smooth_distance(dx0, dy0, dz0, smooth)
                r_1 = 1 / r1
                r_3 = r_1 * r_1 * r_1
                if l >= part_num:
                    m1 = Particles[l, 6]
                    r3 = r_3 * m1
                    A[k, 0] += - dx0 * r3 * (1 + 2 * gamma * m1 * r_1)
                    A[k, 1] += - dy0 * r3 * (1 + 2 * gamma * m1 * r_1)
                    A[k, 2] += - dz0 * r3 * (1 + 2 * gamma * m1 * r_1)
                    A[k, 3] += - (m1 * r_1) * (1 + gamma * m1 * r_1)
                r_3 = r_3 * m
                A[l, 0] += dx0 * r_3 * (1 + 2 * gamma * m * r_1)
                A[l, 1] += dy0 * r_3 * (1 + 2 * gamma * m * r_1)
                A[l, 2] += dz0 * r_3 * (1 + 2 * gamma * m * r_1)
                A[l, 3] += - (m * r_1) * (1 + gamma * m * r_1)
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
