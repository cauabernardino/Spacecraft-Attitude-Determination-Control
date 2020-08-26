from scipy.spatial.transform import Rotation
import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos


### Given constants and other configuration

mrp = np.array([0.1, 0.2, -0.1]).T
mrp_history = [mrp]
w = np.array([np.deg2rad(i) for i in [30, 10, -20]]).T
w_history = [w]

K = 5 # Nm
P = 10 * np.eye(3) #NMs
L = np.array([0.5, -0.3, 0.2]).T #Nm

# Unknown external torque if exists
delta_L = np.array([0.5, -0.3, 0.2]).T #Nm 

I1 = 100
I2 = 75
I3 = 80 # kg m^2

I = np.array([[I1, 0, 0], [0, I2, 0],[0, 0, I3]])


# True for bringing the control to [0, 0, 0] || False for a defined target rate
return_0 = False

###

### Utility functions

def tilde(x):
    x = np.squeeze(x)
    return np.array([[0, -x[2], x[1]],
                     [x[2], 0, -x[0]],
                     [-x[1], x[0], 0]
                     ])


def mrp_to_rotation_matrix(sigma):
    sigma_squared = np.inner(sigma, sigma)
    q0 = (1 - sigma_squared) / (1 + sigma_squared)
    q = [2 * sigma_i / ( 1 + sigma_squared) for sigma_i in sigma]
    q.extend([q0])
    return Rotation.from_quat(q).as_matrix().T

def rotmat_to_mrp(matrix):
    zeta = np.sqrt(np.trace(matrix) + 1)
    constant = 1 / (zeta**2 + 2 * zeta)
    s1 = constant * (matrix[1, 2] - matrix[2, 1])
    s2 = constant * (matrix[2, 0] - matrix[0, 2])
    s3 = constant * (matrix[0, 1] - matrix[1, 0])
    return np.array([s1, s2, s3])

def mrp_shadow(mrp):
    norm = np.linalg.norm(mrp) ** 2
    return np.array([-i / norm for i in mrp])


def target(tval):
    if return_0:
        return np.array([0, 0, 0])
    f = 0.05
    s1 = 0.2 * np.sin(f * tval)
    s2 = 0.3 * np.cos(f * tval)
    s3 = -0.3 * np.sin(f * tval)
    return np.array([s1, s2, s3])

def target_rate(tval):
    if return_0:
        return np.array([0, 0, 0])
    f = 0.05
    s1_dot = 0.2 * f * cos(f * tval)
    s2_dot = -0.3 * f * sin(f * tval)
    s3_dot = -0.3 * f * cos(f * tval)
    sigma_dot = np.array([s1_dot, s2_dot, s3_dot])
    sigma = target(tval)
    A = mrp_dot_matrix(sigma)
    w = 4 * np.dot(np.linalg.inv(A), sigma_dot)

    return w

def target_rate_rate(tval, dt):
    if return_0:
        return np.array([0, 0, 0])
    w1 = target_rate(tval)
    w2 = target_rate(tval - dt)
    return (w1 - w2)/dt

def mrp_dot(mrp, w):
    return 0.25 * np.dot(((1 - np.dot(mrp, mrp)) * np.eye(3) + 2 * tilde(mrp) + 2 * np.outer(mrp,  mrp)), w)

def mrp_dot_matrix(mrp):
    ss = np.dot(mrp, mrp)
    A = np.zeros((3,3))
    A[0, 0] = 1 - ss + 2 * mrp[0] **2
    A[1, 0] = 2*(mrp[1] * mrp[0] + mrp[2])
    A[2, 0] = 2*(mrp[2] * mrp[0] - mrp[1])
    A[0, 1] = 2*(mrp[0] * mrp[1] - mrp[2])
    A[1, 1] = 1 - ss + 2 * mrp[1] ** 2
    A[2, 1] = 2*(mrp[2] * mrp[1] + mrp[0])
    A[0, 2] = 2*(mrp[0] * mrp[2] + mrp[1])
    A[1, 2] = 2*(mrp[1] * mrp[2] - mrp[0])
    A[2, 2] = 1 - ss + 2 * mrp[2] ** 2
    return A

# control u
def control(t, dt, mrp, w):
    sigma_r_n = target(t)
    sigma_b_r = rotmat_to_mrp(np.dot(mrp_to_rotation_matrix(mrp), mrp_to_rotation_matrix(sigma_r_n).T))
    w_r_n = target_rate(t)
    w_r_n_dot = target_rate_rate(t, dt)
    DCM_b_r = mrp_to_rotation_matrix(sigma_b_r)
    w_r_n_body_frame = np.dot(DCM_b_r, w_r_n)
    w_r_n_dot_body_frame = np.dot(DCM_b_r, w_r_n_dot)
    w_b_r = w - w_r_n_body_frame
    """
    EDIT CONTROL HERE!
    """
    # control is
    # -K*mrp_br - P * w_br + I * (w_rn_dot - w_bn X w_rn) + w_bn X Iw_bn
    u = -K * sigma_b_r - np.dot(P, w_b_r) + np.dot(I, w_r_n_dot_body_frame - np.cross(w, w_r_n_body_frame)) + np.cross(w, np.dot(I, w)) - L
    return u

# EOM for w_dot
def wdot(t, dt, mrp, w):
    u = control(t, dt, mrp, w)
    # Iw_dot = -w X Iw + Q

    # Without delta_L
    w_dot = np.dot(np.linalg.inv(I), (-np.cross(w, np.dot(I, w)) + u + L))

    # With delta_L
    # w_dot = np.dot(np.linalg.inv(I), (-np.cross(w, np.dot(I, w)) + u + L + delta_L))
    return w_dot



h = 0.01
time = 100
tvec = np.linspace(0, time, int(time/h + 1))
prev_t = 0
target_histories = [target(0)]
target_rate_history = [target_rate(0)]
error_history = [rotmat_to_mrp(np.dot(mrp_to_rotation_matrix(mrp), mrp_to_rotation_matrix(target(0)).T))]

for ti in tvec[1:]:
    dt = ti - prev_t
    prev_t = ti
    sigma_r_n = target(ti)
    target_histories.append(sigma_r_n)
    sigma_b_r = rotmat_to_mrp(np.dot(mrp_to_rotation_matrix(mrp), mrp_to_rotation_matrix(sigma_r_n).T))
    error_history.append(sigma_b_r)
    w_r_n = target_rate(ti)
    target_rate_history.append(w_r_n)
    w_r_n_dot = target_rate_rate(ti, dt)
    w_b_r = w - w_r_n

    # calculate and apply dots
    mrp = mrp + mrp_dot(mrp, w) * dt
    w = w + wdot(ti, dt, mrp, w) * dt

    if np.dot(mrp, mrp) > 1:
        mrp = mrp_shadow(mrp)

    mrp_history.append(mrp)
    w_history.append(w)
    if ti % 20 == 0:
        print(f"Simulated {ti} seconds")



# Print parameters
moment = 4500 # Seconds * 100

# Print of MRPs
mrp_history = np.array(mrp_history)
mrp_norm = [np.dot(i, i) for i in mrp_history]
print("Norm at ", tvec[moment], "s: ", np.sqrt(mrp_norm[moment]))

w_history = np.array(w_history)
target_histories = np.array(target_histories)
target_rate_history = np.array(target_rate_history)
error_history = np.array(error_history)
# Print of tracking error
error_norm = [np.sqrt(i**2 + j**2 + k**2) for i, j, k in error_history]
print("Norm sigma_br at ", tvec[moment], "s: ", error_norm[moment])

plt.figure(0)
plt.plot(tvec, mrp_history[:, 0], 'g')
plt.plot(tvec, target_histories[:, 0], 'g--')
plt.plot(tvec, mrp_history[:, 1], 'b')
plt.plot(tvec, target_histories[:, 1], 'b--')
plt.plot(tvec, mrp_history[:, 2], 'r')
plt.plot(tvec, target_histories[:, 2], 'r--')
plt.plot(tvec, mrp_norm, 'k')
plt.title('Attitude (sigma) history')
plt.legend(["sigma_1", "sigma_1_target", "sigma_2", "sigma_2_target", "sigma_3", "sigma_3_target", "norm^2"])
plt.grid()
plt.figure(1)
plt.plot(tvec, w_history[:, 0], 'b')
plt.plot(tvec, target_rate_history[:, 0], 'b--')
plt.plot(tvec, w_history[:, 1], 'r')
plt.plot(tvec, target_rate_history[:, 1], 'r--')
plt.plot(tvec, w_history[:, 2], 'g')
plt.plot(tvec, target_rate_history[:, 2], 'g--')
plt.title("Rate (w) history")
plt.legend(["w_1", "w_1_target", "w_2", "w_2_target", "w_3", "w_3_target"])
plt.grid()
plt.show()
