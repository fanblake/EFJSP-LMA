import numpy as np
from FJSP_agv import FJSP_agv,Job,Machine,Agv
import random
import copy
from Data_transfer import Get_Problem
from myplot import save_plot_main_data
import math

def bi_VGM(Pop_size):
    delta=1/Pop_size
    w=[]
    w1=0
    while w1<=1:
        w2=1-w1
        w.append([w1,w2])
        w1+=delta
    return w

# 根据权重向量λ计算T个邻居存入B
# return B->二维list [第一维代表种群中的个体list 第二维代表个体的T个邻居list]
def Neighbor(lambd, T):
    B = []
    for i in range(len(lambd)):
        temp = []
        for j in range(len(lambd)):
            # 计算二维欧氏距离
            distance = np.sqrt((lambd[i][0] - lambd[j][0])**2 + (lambd[i][1] - lambd[j][1])**2)
            temp.append(distance)
        res = np.argsort(temp)  # 下标排序
        B.append(res[:T])  # 取前T个近的邻居加入B
    return B#包含了自己这个向量

def Tchebycheff(value1,value2, z, lambd):
    '''
    :param x: Popi
    :param z: the reference point
    :param lambd: a weight vector
    :return: Tchebycheff objective
    '''
    Gte = []
    Gte.append(np.abs(value1 - z[0]) * lambd[0])
    Gte.append(np.abs(value2 - z[1]) * lambd[0])
    return np.max(Gte)
# 拥挤度选择
def crowding_distance(values1, values2, I):
    distance = [0 for i in range(0, len(I))]
    # sorted1 = sort_by_values(I, values1[:])
    # sorted2 = sort_by_values(I, values2[:])
    sorted1 = sorted(I, key=lambda x: values1[x])  
    sorted2 = sorted(I, key=lambda x: values2[x]) 
    distance[0] = 4444444444444444
    distance[len(I) - 1] = 4444444444444444
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values1[sorted1[k + 1]] - values1[sorted1[k - 1]]) / (
                max(values1) - min(values1))
    for k in range(1, len(I) - 1):
        distance[k] = distance[k] + (values2[sorted2[k + 1]] - values2[sorted2[k - 1]]) / (
                max(values2) - min(values2))#这边我进行了修整，不知道对不对，跟参考的代码编写不同，我感觉我的是对的
    return distance


# Function to sort by values
def sort_by_values(list1, values):
    sorted_list = []
    while (len(sorted_list) != len(list1)):
        if index_of(min(values), values) in list1:
            sorted_list.append(index_of(min(values), values))
        values[index_of(min(values), values)] = math.inf
    return sorted_list


# Function to find index of list
def index_of(a, list):
    for i in range(0, len(list)):
        if list[i] == a:
            return i
    return -1
def fast_non_dominated_sort(values1, values2):
    S = [[] for i in range(0, len(values1))]
    front = [[]]
    n = [0 for i in range(0, len(values1))]
    rank = [0 for i in range(0, len(values1))]

    for p in range(0, len(values1)):
        S[p] = []
        n[p] = 0
        for q in range(0, len(values1)):
            if (values1[p] <= values1[q] and values2[p] <= values2[q]):
                if q not in S[p]:
                    S[p].append(q)
            elif (values1[q] <= values1[p] and values2[q] <= values2[p]):
                n[p] = n[p] + 1
        if n[p] == 0:
            rank[p] = 0
            if p not in front[0]:
                front[0].append(p)

    i = 0
    while (front[i] != []):
        Q = []
        for p in front[i]:
            for q in S[p]:
                n[q] = n[q] - 1
                if (n[q] == 0):
                    rank[q] = i + 1
                    if q not in Q:
                        Q.append(q)
        i = i + 1
        front.append(Q)

    del front[len(front) - 1]
    return front

class MOEAD:
    def __init__(self,Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle,pop_size=200,pc=0.9,pm=0.3,T=5):
        self.fja=FJSP_agv(Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle)
        self.An=An
        self.Pop_size=pop_size
        self.pc=pc
        self.pm=pm
        self.op_num=[len(Pi) for Pi in self.fja.Pt]
        self.Chromo_list=[]
        for i in range(len (self.op_num)):
            self.Chromo_list.extend([i for _ in range(2*self.op_num[i])])
        self.machine_use = [[] for i in range(len(self.op_num))]
        for j in range(len(self.op_num)):
            for k in range(self.op_num[j]):
                mk = []
                for m in range(self.fja.Mn):
                    if self.fja.Pt[j][k][m] != 9999:
                        mk.append(m)
                self.machine_use[j].append(mk)
        self.T=T
        self._lambda = bi_VGM(self.Pop_size)

    def initial_population(self):#随机初始化规则
        Pop1 = []
        Pop2 = []
        Pop3 = []
        # 工件的初始化规则：
        for i in range(self.Pop_size):
            random.shuffle(self.Chromo_list)
            Pop1.append(copy.copy(self.Chromo_list))
            machine_order = []
            for j in range(len(self.machine_use)):
                for k in range(len(self.machine_use[j])):
                    rand_num = random.choice(self.machine_use[j][k])
                    machine_order.append(rand_num)
            Pop2.append(machine_order)
            AGV_order = [random.randint(0, self.An - 1) for i in range(sum(self.op_num))]
            Pop3.append(AGV_order)
        return Pop1,Pop2,Pop3

    def jbchrom_reset(self,pop1,pop3):#调整第一层染色体的顺序，防止AGV的载重超载
        # for Pi in range(len(pop1)):
        c = [0 for i in range(self.fja.Jn)]  # 用来记录装卸及工序号
        c_ = [0 for i in range(self.fja.Jn)]  # 用来记录工序的数量，只有卸载之后才能+1
        a = [0 for i in range(self.fja.An)]  # 用来记录AGV的载重
        a1 = [[] for i in range(self.fja.An)]  # 用来记录AGV的加工工件的路径
        a2 = [[] for i in range(self.fja.An)]  # 用来记录基因的位置
        b = 0  # 需要删除基因的位置
        for Pj in range(len(pop1)):
            j=pop1[Pj]#工件号
            if j==0:
                gx_loca=c_[j]#工序号
            else:#不是第一个工件
                gx_loca=sum(self.op_num[:j])+c_[j]#工序的位置
            agv=pop3[gx_loca]#agv号
            if c[j]%2==0: #表示该工序是装载，需要判断载重,且为第一道工序
                if a[agv]+1>AGV_weight:#超过了载重，需要更新基因的位置以及机器的位置
                    #找到上一个加工的位置的基因 然后调整到目前这个基因的位置前
                    #先更换位置，然后卸载工序
                    j=a1[agv][-1]#上一个工件号
                    for i in range(a2[agv][-1]+1,len(pop1)):
                        if pop1[i]==j:
                            b=i
                            # print("调整染色体")
                            break
                    pop1.pop(b)#删除该位置的索引的值
                    pop1.insert(Pj,j)
                    c_[j] += 1
                    a[agv] -= 1
                else:
                    a[agv]+=1
            else:#表示该工序是卸载,需要更新AGV的载重
                c_[j]+=1
                a[agv]-=1
            c[j] += 1
            a1[agv].append(j)
            a2[agv].append(Pj)
        return pop1
    def cross(self,pop1,pop2,pop3,j1,j2):
        new_pop1 = [-1 for i in range(len(pop1[0]))]
        new_pop2 = [-1 for i in range(len(pop2[0]))]
        new_pop3 = [-1 for i in range(len(pop3[0]))]
        gj = [i for i in range(self.fja.Jn)]
        # if np.random.rand() < self.pc:
        random.shuffle(gj)  # MPX交叉
        # 工序交叉
        split_point = random.randint(1, len(gj) - 1)
        x1 = gj[:split_point]
        x2 = gj[split_point:]
        for j in range(len(pop1[j1])):
            if pop1[j1][j] in x1:
                new_pop1[j] = pop1[j1][j]
        site = [i for i, var in enumerate(new_pop1) if var == -1]
        e = 0
        for j in range(len(pop1[j1])):
            if pop1[j2][j] in x2:
                new_pop1[site[e]] = pop1[j2][j]
                e = e + 1
        # 机器与AGV的交叉
        random_list = [random.choice([0, 1]) for _ in range(sum(self.op_num))]
        for j in range(len(random_list)):
            if random_list[j] == 1:
                new_pop2[j] = pop2[j1][j]
                new_pop3[j] = pop3[j1][j]
            else:
                new_pop2[j] = pop2[j2][j]
                new_pop3[j] = pop3[j2][j]
        return new_pop1,new_pop2,new_pop3


    def mutation(self,pop1,pop2,pop3):
        new_pop1 = copy.copy(pop1)
        new_pop2 = copy.copy(pop2)
        new_pop3 = copy.copy(pop3)
        if np.random.rand() <self.pm:
            random_positions=random.sample(range(len(pop1)),2)
            b=new_pop1[random_positions[0]]
            new_pop1[random_positions[0]]=new_pop1[random_positions[1]]
            new_pop1[random_positions[1]]=b
            random_positions1=random.sample(range(len(pop3)),3)
            for j in random_positions1:
                # new_pop2[j]=random.choice[i for i in range(self.Machine) if i=new_pop2[j] and ]#不更换机器了，比较复杂
                new_pop3[j]=random.choice([s for s in range(self.An) if s!=new_pop3[j]])
        return new_pop1,new_pop2,new_pop3

    def Dominate(self, Pop1, Pop2):
        '''
        :param Pop1:
        :param Pop2:
        :return: If Pop1 dominate Pop2, return True
        '''
        func_value1=[]
        func_value2=[]

        fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk,
                       epik, ea_load, ea_unload, ea_idle)
        self.Jobs, self.Machines, self.AGVs = fja.decode(Pop1[0], Pop1[1], Pop1[2])
        v1, v2 = fja.total_value()
        func_value1.append(v1)
        func_value2.append(v2)
        fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk,
                       epik, ea_load, ea_unload, ea_idle)
        self.Jobs, self.Machines, self.AGVs = fja.decode(Pop2[0], Pop2[1], Pop2[2])
        v1, v2 = fja.total_value()
        func_value1.append(v1)
        func_value2.append(v2)
        if (func_value1[0]<func_value1[1] and func_value2[0]<func_value2[1]) or \
                (func_value1[0] <= func_value1[1] and func_value2[0] < func_value2[1]) or \
                (func_value1[0] < func_value1[1] and func_value2[0] <= func_value2[1]):
            return True
        else:
            return False




    def MOEAD_main(self,pop1,pop2,pop3):
        # 获得邻域向量
        B = Neighbor(MO._lambda, MO.T)  # work out the T closest weight vectors to each weight vector
        # 存储非支配解
        EP = []  # EP is used to store non-dominated solutions found during the search
        z=[]
        func_value1 = []
        func_value2 = []
        for i in range(len(pop1)):
            fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk,
                           epik, ea_load, ea_unload, ea_idle)
            self.Jobs, self.Machines, self.AGVs = fja.decode(pop1[i], pop2[i], pop3[i])
            v1, v2 = fja.total_value()
            func_value1.append(v1)
            func_value2.append(v2)
        z.append(min(func_value1))
        z.append(min(func_value2))#参考点
        pareto_front_solution = []
        func1_best = []
        func2_best = []
        for g in range(GEN):#迭代次数
            print('这是算法的第几次迭代*****************************************:', g)
            for i in range(len(pop1)):
                # Randomly select two indexes k,l from B(i)
                j = random.randint(0, self.T - 1)
                k = random.randint(0, self.T - 1)
                new_pop1,new_pop2,new_pop3=MO.cross(pop1,pop2,pop3,B[i][j], B[i][k])
                new_pop1 = MO.jbchrom_reset(new_pop1, new_pop3)  # 修整工序染色体
                # new_pop1, new_pop2, new_pop3 = MO.mutation(new_pop1,new_pop2,new_pop3)#本身该产生两个个体，我只采用了其中一个个体
                # new_pop1 = MO.jbchrom_reset(new_pop1, new_pop3)  # 修整工序染色体
                fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk,
                               epik, ea_load, ea_unload, ea_idle)
                self.Jobs, self.Machines, self.AGVs = fja.decode(new_pop1, new_pop2, new_pop3)
                v1, v2 = fja.total_value()
                if z[0]>v1:
                    z[0]=v1
                if z[1]>v2:
                    z[1]=v2  #更新参考点
                # update of Neighboring solutions
                for bi in range(len(B[i])):
                    Ta=Tchebycheff(func_value1[B[i][bi]],func_value2[B[i][bi]], z, self._lambda[B[i][bi]])
                    Tb = Tchebycheff(v1,v2, z, self._lambda[B[i][bi]])
                    if Tb < Ta:
                        pop1[B[i][j]],pop2[B[i][j]],pop3[B[i][j]] = new_pop1,new_pop2,new_pop3
                        func_value1[B[i][j]]=v1
                        func_value2[B[i][j]]=v2
                y=[]
                y.append(new_pop1)
                y.append(new_pop2)
                y.append(new_pop3)
                # Update of EP
                if EP == []:
                    EP.append(y)
                else:
                    dominateY = False  # 是否有支配Y的解
                    _remove = []  # Remove from EP all the vectors dominated by y
                    for ei in range(len(EP)):
                        if MO.Dominate(y, EP[ei]):
                            _remove.append(EP[ei])
                        elif MO.Dominate(EP[ei], y):
                            dominateY = True
                            break
                    # add y to EP if no vectors in EP dominated y
                    if not dominateY:
                        EP.append(y)
                        for j in range(len(_remove)):
                            EP.remove(_remove[j])
                # 1、从EP中移除被y支配的向量，2、如果没有向量支配y，将y加入EP
            print(len(EP))
            f_value1 = []
            f_value2 = []
            for i in range(len(EP)):
                fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight,
                               epk,
                               epik, ea_load, ea_unload, ea_idle)
                self.Jobs, self.Machines, self.AGVs = fja.decode(EP[i][0], EP[i][1], EP[i][2])
                v1, v2 = fja.total_value()
                f_value1.append(v1)
                f_value2.append(v2)
            objectives = []
            for j in range(len(EP)):
                object = []
                object.append(f_value1[j])
                object.append(f_value2[j])
                objectives.append(object)
            pareto_front_solution.append(objectives)
            ob = np.asarray(objectives)
            # func1_best.append(sum(ob[:,0])//len(ob))
            # func2_best.append(sum(ob[:, 1])//len(ob))
            func1_best.append(min(ob[:, 0]))
            func2_best.append(min(ob[:, 1]))
            print(pareto_front_solution)
            print(func1_best)
        return pareto_front_solution,func1_best,func2_best#所有迭代完成后的精英解



if __name__=="__main__":
    # files = [ 11, 13, 15]
    # for file_index in range(7,16):
    #     print(f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')
    #     P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk, epik, ea_load,ea_unload, ea_idle= Get_Problem(
    #         f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')
    GEN = 400
    ps1 = 100
    pc1 = 0.8
    pm1 = 0.2

    Machine_number = 10
    Job_number = 7
    AGV_number = 3
    AGV_weight = 2
    # epk = [random.randint(10, 18) for i in range(Machine_number)]
    epk = [10, 8, 12, 10, 7, 13, 11, 8, 10, 11]
    # epik = [random.randint(1, 4) for i in range(Machine_number)]
    epik = [2, 1, 2, 1, 2, 1, 1, 2, 3, 2]
    ea_unload = 3.5
    ea_load = 0.5
    ea_idle = 0.5
    P_time = [[
        [9999, 9999, 6, 4, 9999, 9999, 9999, 9999, 9999, 9999],
        [2, 3, 9999, 9999, 4, 9999, 9999, 9999, 9999, 9999],
        [4, 3, 3, 6, 9999, 4, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 4, 3, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 3, 5]
    ], [
        [9999, 9999, 3, 4, 9999, 9999, 9999, 9999, 9999, 9999],
        [2, 2, 4, 3, 9999, 9999, 9999, 9999, 9999, 9999],
        [3, 2, 9999, 9999, 4, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 4, 6, 9999, 9999, 9999, 9999, 9999, 9999],
        [3, 3, 9999, 9999, 9999, 9999, 9999, 9999, 2, 4],
        [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 7, 8]

    ], [
        [5, 7, 9999, 9999, 9999, 9999, 9999, 9999, 9, 6],
        [3, 5, 9999, 9999, 3, 9999, 9999, 9999, 9999, 9999],
        [7, 5, 9999, 9999, 9999, 4, 9999, 9999, 9999, 9999],
        [3, 2, 2, 4, 9999, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 6, 9, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 4, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 4, 7, 9999, 9999]
    ], [
        [4, 4, 3, 6, 9999, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9, 7],
        [2, 1, 9999, 9999, 3, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 7, 6, 7, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 4, 9999, 9999, 9999, 9999],
        [3, 3, 9999, 9999, 9999, 9999, 9999, 9999, 4, 5]
    ], [
        [9999, 9999, 5, 6, 5, 9999, 9999, 9999, 9999, 9999],
        [3, 2, 9999, 9999, 2, 9999, 9999, 9999, 9999, 9999],
        [6, 6, 5, 8, 5, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 2, 3],
        [9999, 9999, 9999, 9999, 9999, 9999, 4, 3, 9999, 9999],
        [3, 2, 9999, 9999, 2, 9999, 9999, 9999, 9999, 9999]

    ], [
        [3, 4, 9999, 9999, 9999, 9999, 9999, 9999, 5, 5],
        [9999, 9999, 9999, 9999, 9999, 9999, 4, 3, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 3, 9999, 9999, 9999, 9999],
        [9999, 9999, 3, 4, 3, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 4, 6],
        [3, 3, 2, 4, 9999, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 3, 5, 9999, 9999]
    ], [
        [3, 3, 9999, 2, 9999, 9999, 9999, 9999, 4, 5],
        [3, 2, 9999, 9999, 3, 9999, 9999, 9999, 9999, 9999],
        [1, 9999, 2, 1, 9999, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 4, 4, 3, 9999, 9999, 9999, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 3, 4, 9999, 9999],
        [9999, 9999, 9999, 9999, 9999, 9999, 5, 4, 9999, 9999],
        [3, 4, 9999, 4, 9999, 9999, 9999, 9999, 3, 5]
    ]]

    AGV_Trans = [[0, 1.5, 3.9, 3.1, 3.3, 4.1, 4.6, 3.3, 3.7, 3.0, 3.3],
                 [1.5, 0, 2.4, 1.6, 1.8, 2.4, 3.1, 1.8, 2.2, 1.5, 1.8],
                 [3.9, 2.4, 0, 1.7, 1.5, 2.3, 2.2, 1.2, 1.0, 2.5, 2.2],
                 [3.1, 1.6, 1.7, 0, 1.0, 2.2, 2.8, 1.9, 1.2, 1.4, 1.2],
                 [3.3, 1.8, 1.5, 1.0, 0, 2.7, 2.4, 1.8, 1.6, 2.0, 1.8],
                 [4.1, 2.4, 2.3, 2.2, 2.7, 0, 1.8, 2.8, 2.6, 1.2, 3.0],
                 [4.6, 3.1, 2.2, 2.8, 2.4, 1.8, 0, 2.0, 2.3, 2.8, 2.6],
                 [3.3, 1.8, 1.2, 1.9, 1.8, 2.8, 2.0, 0, 1.0, 1.2, 1.4],
                 [3.7, 2.2, 1.0, 1.2, 1.6, 2.6, 2.3, 1.0, 0, 1.5, 1.7],
                 [3.0, 1.5, 2.5, 1.4, 2.0, 1.2, 2.3, 1.2, 1.5, 0, 1.1],
                 [3.3, 1.8, 2.2, 1.2, 1.8, 3.0, 2.6, 1.4, 1.7, 1.1, 0]
                 ]

    MO=MOEAD(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk, epik, ea_load,
            ea_unload, ea_idle)
    pop1, pop2, pop3 =MO.initial_population()  # 初始化
    for p in range(len(pop1)):
        pop1[p] = MO.jbchrom_reset(pop1[p], pop3[p])  # 修整工序染色体
    pareto_front_solution, func1_best, func2_best = MO.MOEAD_main(pop1, pop2, pop3)
    save_plot_main_data("pic_data/MOEAD/data/",  pareto_front_solution,func1_best,func2_best,GEN)
    #这个也是生成两个个体