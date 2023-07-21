import numpy as np
import sqtdiat.qops as sq
import matplotlib.pyplot as plt
from qiskit.visualization import array_to_latex
import seaborn as sns

def _make_unitary(gate, pos, num):
    """
    Creates a unitary qutrit gate with the position as a user input.

    Parameters:
    - gate: The unitary gate to be placed at the specified position. It should be a 3x3 NumPy array representing the gate's transformation.
    - pos: The position where the gate should be placed within the larger unitary gate. The position is zero-indexed.
    - num: The total number of qutrits in the system. It represents the size of the larger unitary gate to be created.

    Returns:
    - ini_gate: The resulting unitary gate, represented as a NumPy array, where the specified gate is placed at the specified position within a larger unitary gate.
    """

    if pos > num - 1:
            raise ValueError(f"{pos} not within {num}")
    else:
        if pos == 0:
            ini_gate = gate
            for i in range(num - 1):
                ini_gate = np.kron(ini_gate, np.eye(3))
        else:
            ini_gate = np.eye(3)
            for i in range(num - 1):
                if pos == i + 1:
                    ini_gate = np.kron(ini_gate, gate)
                else:
                    ini_gate = np.kron(ini_gate, np.eye(3))
        return ini_gate
    
def d2t(num, rng):
    """
    Converts a decimal number to a ternary (base-3) representation with a specified range.

    Parameters:
    - num: The decimal number to be converted to ternary representation.
    - rng: The desired length of the ternary representation. If the converted ternary representation is shorter than the specified range, leading zeros will be added.

    Returns:
    - gt: The ternary representation of the decimal number as a string, with leading zeros if necessary.
    """

    gt = ''
    while num > 0:
        gt = gt + str(num % 3)
        num //= 3
    while len(gt) < rng:
        gt = gt + '0'
    return gt[::-1]

def d2tarr(num, rng):
    """
    Converts a decimal number to a ternary array representation.

    Parameters:
    - num: The decimal number to be converted.
    - rng: The desired length of the ternary array representation.

    Returns:
    - mm: A list representing the ternary array obtained from the conversion of the decimal number. The list will have a length equal to the specified range (rng).
    """
    mm = []
    while num > 0:
        mm.append(num % 3)
        num //= 3
    while len(mm) < rng:
        mm.append(0)
    return mm[::-1]

def t2d(str_):
    num = 0
    str_ = str_[::-1]
    for i in range(len(str_)):
        num += int(str_[i]) * 3**i
    return num

def cx_replace(ctr, tgt, num, rng):
    ter = d2tarr(num, rng)
    ter[tgt] = (ter[tgt] - ter[ctr]) % 3
    return t2d(ter)

def _make_ctrl_uni(num, gate_, ctrl, trgt, pos_num):
    den = np.eye(3**num) + 1j * np.zeros((3**num, 3**num))
    if len(pos_num) == 1:
        for i in range(3**num):
            t_val = d2t(i, num)
            if t_val[ctrl] == str(pos_num[0]):
                t_val_arr = [int(x) for x in t_val]
                t_val_arr[trgt] = 0
                for j in range(3):
                    den[t2d(t_val), t2d(t_val_arr)] = gate_[t2d(str(t_val_arr[trgt])), t2d(str(t_val[trgt]))]
        #             print(t_val, t_val_arr)
                    t_val_arr[trgt] += 1
    elif len(pos_num) == 2:
        gate_2 = gate_ @ gate_
        for position in pos_num:
            for i in range(3**num):
                t_val = d2t(i, num)
                if t_val[ctrl] == '1':
                    t_val_arr = [int(x) for x in t_val]
                    t_val_arr[trgt] = 0
                    for j in range(3):
                        den[t2d(t_val), t2d(t_val_arr)] = gate_[t2d(str(t_val_arr[trgt])), t2d(str(t_val[trgt]))]
            #             print(t_val, t_val_arr)
                        t_val_arr[trgt] += 1
                elif t_val[ctrl] == '2':
                    t_val_arr = [int(x) for x in t_val]
                    t_val_arr[trgt] = 0
                    for j in range(3):
                        den[t2d(t_val), t2d(t_val_arr)] = gate_2[t2d(str(t_val_arr[trgt])), t2d(str(t_val[trgt]))]
#                         print(t_val, t_val_arr)
                        t_val_arr[trgt] += 1
    return den

class Qutrit:
    def __init__(self, num):
        """
        Initializes a Qutrit object.

        Parameters:
        - num: The number of qutrits in the system.
        """
        self.num = num
        state = np.array([1, 0, 0])
        for i in range(self.num - 1):
            state = np.kron(state, np.array([1, 0, 0]))
        self.state = state
        self.circ_store = [[] for i in range(num)]
        
    def initialize(self, arr, pos):
        """
        Initializes the qutrit state with a specified array at the given position.

        Parameters:
        - arr: An array representing the desired state.
        - pos: The position where the state should be initialized within the qutrit(s).
        """
        if pos == 0:
            ini_state = np.array(arr)
            for i in range(self.num - 1):
                ini_state = np.kron(ini_state, np.array([1, 0, 0]))
        else:
            ini_state = np.array([1, 0, 0])
            for i in range(self.num - 1):
                if pos == i + 1:
                    ini_state = np.kron(ini_state, np.array(arr))
                else:
                    ini_state = np.kron(ini_state, np.array([1, 0, 0]))
        self.state = ini_state
        
    def _draw_ctr(self, gate_name, ctr, trgt):
        max_ = []
        for x in self.circ_store: max_.append(len(x)) 
        for i in range(self.num): 
            while len(self.circ_store[i]) <= np.max(max_): 
                self.circ_store[i].append('---')
        for i in range(self.num):
            if i == ctr:
                self.circ_store[i].append(' * ')
            elif i == trgt:
                self.circ_store[i].append(gate_name)
            elif ctr < trgt and ctr < i < trgt:
                self.circ_store[i].append('-|-')
            elif ctr > trgt and trgt < i < ctr:
                self.circ_store[i].append('-|-')
            else:
                self.circ_store[i].append('---')
            
    def plot_density(self, p_q = [],  axis = False, x = 15, y = 6, an = True):
        """
        Plots the density matrix of the qutrit.

        Parameters:
        - axis: Boolean flag to show/hide the axis. Default is False.
        - x: The width of the plot. Default is 15.
        - y: The height of the plot. Default is 6.
        - an: Boolean flag to show/hide the annotations. Default is True.
        """
        if len(p_q) == 0:
            den = self.density()
        else:
            den = Qutrit.partial_trace(self, p_q)
        plt.figure(figsize = (x, y))
        plt.subplot(1, 2, 1)
        sns.heatmap(np.round(np.real(den), 5), annot = an)
        plt.title("Real part")
        if axis == False:
            plt.axis('off')
        plt.subplot(1, 2, 2)
        sns.heatmap(np.round(np.imag(den), 5), annot = an)
        plt.title('Imaginary part')
        if axis == False:
            plt.axis('off')
        plt.show()
            
    def show_statevector(self):
        """
        Returns the LaTeX representation of the qutrit's state vector.
        """
        return array_to_latex(self.state)
    
    def density(self):
        return np.outer(self.state, np.conj(self.state))
    
    def partial_trace(self, pttx):
        """
        Made originally by Umesh Chandra Joshi, M.Tech(Quantum Computing) 2022-24 Batch, D.I.A.T, India, distributed under MIT License.
        https://github.com/UMESHJOSHI106, Contact : umesh3210joshi@gmail.com
        Returns the partial trace of the qutrit's density matrix.

        Parameters:
        - pttx: A list of indices specifying the qutrits to trace out.
        
        """
        rho = self.density()
        tot = self.num
        kk = [3**i for i in pttx]
        kk = kk[::-1]
        jj = [3**i for i in range(tot) if i not in pttx]
        jj = jj[::-1]
        
        def binfromdec(number, rng):
            arr = []
            while number > 0:
                arr.append(str(number % 3))
                number //= 3
            while len(arr) < rng:
                arr.append(str(0))
            return "".join(arr[::-1])


        gg = []

        for i in range(3**(tot - len(pttx))):
            str1 = binfromdec(i, tot - len(pttx))
            sum = 0
            for j in range(len(str1)):
                sum += int(str1[j]) * jj[j]
            gg.append(sum)

        hh = []
        for i in range(3**(len(pttx))):
            str1 = binfromdec(i, len(pttx))
            sum = 0
            for j in range(len(str1)):
                sum += int(str1[j]) * kk[j]
            hh.append(sum)
        hh = np.array(hh)

        ans = np.zeros((3**(tot - len(pttx)), 3**(tot - len(pttx)))) + 1j * np.zeros((3**(tot - len(pttx)), 3**(tot - len(pttx))))
        for i in range(len(hh)):
            ll = hh[i] + gg
            for j in range(ans.shape[0]):
                for k in range(ans.shape[0]):
                    ans[j][k] += rho[ll[j]][ll[k]]
        return ans
    
    def get_probs(self, x_l = 20, p_q = []):
        """
        Computes and plots the probabilities of the qutrit's states.

        Parameters:
        - x_l: The width of the plot. Default is 20.
        - p_q: A list of indices specifying the qutrits to consider for probability computation. Default is an empty list.
        """
        plt.rcParams.update({'font.size': 22})
        keys = []
        vals = []
        if len(p_q) == 0:
            den = self.density()
        else:
            den = Qutrit.partial_trace(self, p_q)
        for i in range(den.shape[0]):
            a = np.real(den[i, i])
            if a > 1e-4:
                vals.append(a)
                keys.append(d2t(i, self.num - len(p_q)))
        plt.figure(figsize = (x_l, 3 * self.num))
        plt.bar(keys, vals)
        plt.ylim(0, max(vals) + 0.1 * max(vals))
        plt.title("Probabilities")
        for i, height in enumerate(vals):
            if height > 1e-4:
                plt.text(i, height, str(np.round(height, 3)), ha='center', va='bottom')
        plt.show()
        
    def draw(self):
        for x in self.circ_store:
            for y in x:
                print(y, end = '')
                print('-', end = '')
            print()
        
    def X01(self, pos):
        x01 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        self.state = _make_unitary(x01, pos, self.num) @ self.state
        self.circ_store[pos].append('X01')
        
    def X12(self, pos):
        x12 = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
        self.state = _make_unitary(x12, pos, self.num) @ self.state
        self.circ_store[pos].append('X12')
    
    def X02(self, pos):
        x02 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        self.state = _make_unitary(x02, pos, self.num) @ self.state
        self.circ_store[pos].append('X02')
    
    def XM1(self, pos):
        xm1 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        self.state = _make_unitary(xm1, pos, self.num) @ self.state
        self.circ_store[pos].append('X-1')
    
    def XP1(self, pos):
        xp1 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        self.state = _make_unitary(xp1, pos, self.num) @ self.state
        self.circ_store[pos].append('X+1')
        
    def ZP1(self, pos):
        om = np.exp(2 * np.pi * 1j / 3)
        zp1 = np.array([[1, 0, 0], [0, om, 0], [0, 0, om**2]])
        self.state = _make_unitary(zp1, pos, self.num) @ self.state
        self.circ_store[pos].append('Z+1')
        
    def ZM1(self, pos):
        om = np.exp(2 * np.pi * 1j / 3)
        zm1 = np.array([[1, 0, 0], [0, om**2, 0], [0, 0, om]])
        self.state = _make_unitary(zm1, pos, self.num) @ self.state
        self.circ_store[pos].append('Z-1')
            
    def H(self, pos):
        om = np.exp(2 * np.pi * 1j / 3)
        h = -1j / np.sqrt(3) * np.array([[1, 1, 1], [1, om, om**2], [1, om**2, om]])
        self.state = _make_unitary(h, pos, self.num) @ self.state
        self.circ_store[pos].append(' H ')
        
    def H_DAG(self, pos):
        om = np.exp(2 * np.pi * 1j / 3)
        h = -1j / np.sqrt(3) * np.array([[1, 1, 1], [1, om, om**2], [1, om**2, om]])
        self.state = _make_unitary(np.conj(h).T, pos, self.num) @ self.state
        self.circ_store[pos].append('H_d')
        
    def S(self, pos):
        """
        Applies the S gate (Phase gate) to the qutrit at the specified position.

        Parameters:
        - pos: The position of the qutrit to apply the gate.
        """
        s = np.diag([1, 1, np.exp(2j * np.pi / 3)])
        self.state = _make_unitary(s, pos, self.num) @ self.state
        self.circ_store[pos].append(' S ')
        
    def Z0(self, pos):
        z0 = np.diag([1, 1, -1])
        self.state = _make_unitary(z0, pos, self.num) @ self.state
        self.circ_store[pos].append('Z_0')
    
    def Z1(self, pos):
        z1 = np.diag([1, np.exp(2j * np.pi / 3), np.exp(-2j * np.pi / 3)])
        self.state = _make_unitary(z1, pos, self.num) @ self.state
        self.circ_store[pos].append('Z_1')
    
    def Z2(self, pos):
        z2 = np.diag([1, np.exp(-2j * np.pi / 3), np.exp(2j * np.pi / 3)])
        self.state = _make_unitary(z2, pos, self.num) @ self.state
        self.circ_store[pos].append('Z_2')
        
    def CX(self, ctr, trgt, show_gate = False):
        if ctr == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctr < self.num and trgt < self.num:
            cx_ = np.eye(3**self.num)
            for i in range(3**self.num):
                j = cx_replace(ctr, trgt, i, self.num)
                cx_[i, i], cx_[i, j] = cx_[i, j], cx_[i, i]
            if show_gate == True:
                return cx_
            self.state = cx_ @ self.state
            Qutrit._draw_ctr(self, '-X-', ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
            
    def CX_DAG(self, ctr, trgt, show_gate = False):
        if ctr == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctr < self.num and trgt < self.num:
            cx_ = np.eye(3**self.num)
            for i in range(3**self.num):
                j = cx_replace(ctr, trgt, i, self.num)
                cx_[i, i], cx_[i, j] = cx_[i, j], cx_[i, i]
            if show_gate == True:
                return np.linalg.pinv(cx_)
            self.state = np.linalg.pinv(cx_) @ self.state
            Qutrit._draw_ctr(self, 'X_d', ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
    
#     def CU_2(self, gate_, ctrl, trgt):
#         if ctrl == trgt:
#             raise ValueError("Control cannot be same as Target")
#         if ctrl < self.num and trgt < self.num:
#             num = self.num
#             den = np.eye(3**num) + 1j * np.zeros((3**num, 3**num))
#             for i in range(3**num):
#                 t_val = d2t(i, num)
#                 if t_val[ctrl] == '2':
#                     t_val_arr = [int(x) for x in t_val]
#                     t_val_arr[trgt] = 0
#                     for j in range(3):
#                         den[t2d(t_val), t2d(t_val_arr)] = gate_[t2d(str(t_val_arr[trgt])), t2d(str(t_val[trgt]))]
#             #             print(t_val, t_val_arr)
#                         t_val_arr[trgt] += 1
#             self.state = den @ self.state
#         else:
#             print(self.num, ctrl, trgt)
#             raise ValueError("Control and Target should be less than the number of qutrits.")
            
    
    def CU_2(self, gate_, ctrl, trgt, show_gate = False):
        if ctrl == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctrl < self.num and trgt < self.num:
            den = _make_ctrl_uni(self.num, gate_, ctrl, trgt, ['2'])
            if show_gate == True:
                return den
            self.state = den @ self.state
            Qutrit._draw_ctr(self, 'CU2', ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
        
    def CU_2_DAG(self, gate_, ctrl, trgt, show_gate = False):
        if ctrl == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctrl < self.num and trgt < self.num:
            den = _make_ctrl_uni(self.num, gate_, ctrl, trgt, ['2'])
            if show_gate == True:
                return np.linalg.pinv(den)
            self.state = np.linalg.pinv(den) @ self.state
            Qutrit._draw_ctr(self, 'C2d', ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
            
    def CU_1(self, gate_, ctrl, trgt, show_gate = False):
        if ctrl == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctrl < self.num and trgt < self.num:
            den = _make_ctrl_uni(self.num, gate_, ctrl, trgt, ['1'])
            if show_gate == True:
                return den
            self.state = den @ self.state
            Qutrit._draw_ctr(self, 'CU1', ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
            
    def CU_1_DAG(self, gate_, ctrl, trgt, show_gate = False, name = ''):
        if ctrl == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctrl < self.num and trgt < self.num:
            den = _make_ctrl_uni(self.num, gate_, ctrl, trgt, ['1'])
            if show_gate == True:
                return np.linalg.pinv(den)
            self.state = np.linalg.pinv(den) @ self.state
            if name == '':
                nm = 'C1d'
            else: name = nm
            Qutrit._draw_ctr(self, nm, ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
            
    def CU(self, gate_, ctrl, trgt, show_gate = False):
        if ctrl == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctrl < self.num and trgt < self.num:
            den = _make_ctrl_uni(self.num, gate_, ctrl, trgt, ['1', '2'])
            if show_gate == True:
                return den
            self.state = den @ self.state
            Qutrit._draw_ctr(self, 'C_U', ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
            
    def CU_DAG(self, gate_, ctrl, trgt, show_gate = False):
        if ctrl == trgt:
            raise ValueError("Control cannot be same as Target")
        if ctrl < self.num and trgt < self.num:
            den = _make_ctrl_uni(self.num, gate_, ctrl, trgt, ['1', '2'])
            if show_gate == True:
                return np.linalg.pinv(den)
            self.state = np.linalg.pinv(den) @ self.state
            Qutrit._draw_ctr(self, 'CUd', ctr, trgt)
        else:
            raise ValueError("Control and Target should be less than the number of qutrits.")
            
#     def CU_DAG(self, gate_, ctrl, trgt, show_gate = False):
#         if ctrl == trgt:
#             raise ValueError("Control cannot be same as Target")
#         if ctrl < self.num and trgt < self.num:
#             num = self.num
#             den = np.eye(3**num) + 1j * np.zeros((3**num, 3**num))
#             gate_2 = gate_ @ gate_
#             for i in range(3**num):
#                 t_val = d2t(i, num)
#                 if t_val[ctrl] == '1':
#                     t_val_arr = [int(x) for x in t_val]
#                     t_val_arr[trgt] = 0
#                     for j in range(3):
#                         den[t2d(t_val), t2d(t_val_arr)] = gate_[t2d(str(t_val_arr[trgt])), t2d(str(t_val[trgt]))]
#             #             print(t_val, t_val_arr)
#                         t_val_arr[trgt] += 1
#                 elif t_val[ctrl] == '2':
#                     t_val_arr = [int(x) for x in t_val]
#                     t_val_arr[trgt] = 0
#                     for j in range(3):
#                         den[t2d(t_val), t2d(t_val_arr)] = gate_2[t2d(str(t_val_arr[trgt])), t2d(str(t_val[trgt]))]
# #                         print(t_val, t_val_arr)
#                         t_val_arr[trgt] += 1
#             if show_gate == True:
#                 return np.linalg.pinv(den)
#             self.state = np.linalg.pinv(den) @ self.state
#         else:
#             raise ValueError("Control and Target should be less than the number of qutrits.")
