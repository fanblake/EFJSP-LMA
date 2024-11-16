#初始化种群
#非支配排序和拥挤度计算
#选择，交叉，变异，锦标赛选择
#种群合并（2N)
#非支配排序与拥挤度计算
#生成新的种群（N）
from FJSP_MlAGV import MA
from Data_transfer import Get_Problem
from FJSP_agv import FJSP_agv,Job,Machine,Agv
import random
import copy
import math
import numpy as np
from myplot import save_plot_main_data

# 快速非支配排序
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


# 拥挤度选择
def crowding_distance(values1, values2, I):
    distance = [0 for i in range(0, len(I))]
    # sorted1 = sort_by_values(I, values1[:])
    # sorted2 = sort_by_values(I, values2[:])
    sorted1 = sorted(I, key=lambda x: values1[x])  
    sorted2 = sorted(I, key=lambda x: values2[x]) 
    distance[0] = 4444444444444444
    distance[len(I) - 1] = 4444444444444444
    if(max(values1) == min(values1)):
        print(values1)
        print(values2)
        print("eeeeeeeeeeeeeeeeeeeeeeeee")
        pass
    else:
        for k in range(1, len(I) - 1):
            distance[k] = distance[k] + (values1[sorted1[k + 1]] - values1[sorted1[k - 1]]) / (
                    max(values1) - min(values1))
        for k in range(1, len(I) - 1):
            distance[k] = distance[k] + (values2[sorted2[k + 1]] - values2[sorted2[k - 1]]) / (
                    max(values2) - min(values2))  # 这边我进行了修整，不知道对不对，跟参考的代码编写不同，我感觉我的是对的
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
class NSGAII:
    def __init__(self,Pt,Jn,Mn,At,An,Aw,epk,epik,ea_load,ea_unload,ea_idle,pop_size=200,pc=0.8,pm=0.2,N_elite=10):
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

    #锦标赛选择
    def Mating(self, pop_num, rank,crowd):#锦标赛选择
        Matingpool=[]
        num_range = list(range(0, pop_num))
        selected_pairs = []  # 用于存储已经被选择的个体
        for b in range(pop_num//2):
            #随机选择两个个体，进行比较，好的个体留下,避免重复选择两个个体
            select_numbers=random.sample(num_range,2)
            selected_pairs.append(select_numbers)
            for num in select_numbers:
                num_range.remove(num)
            p=select_numbers[0]
            q=select_numbers[1]
            if rank[p] < rank[q]:
                best = p
            elif rank[p] > rank[q]:
                best = q
            elif crowd[p] > crowd[q]:
                best = p
            else:
                best = q
            Matingpool.append(best)
        return Matingpool
    def jbchrom_reset(self,pop1,pop3):#调整第一层染色体的顺序，防止AGV的载重超载
        for Pi in range(len(pop1)):
            c = [0 for i in range(self.fja.Jn)]  # 用来记录装卸及工序号
            c_ = [0 for i in range(self.fja.Jn)]  # 用来记录工序的数量，只有卸载之后才能+1
            a = [0 for i in range(self.fja.An)]  # 用来记录AGV的载重
            a1 = [[] for i in range(self.fja.An)]  # 用来记录AGV的加工工件的路径
            a2 = [[] for i in range(self.fja.An)]  # 用来记录基因的位置
            b = 0  # 需要删除基因的位置
            for Pj in range(len(pop1[Pi])):
                j=pop1[Pi][Pj]#工件号
                if j==0:
                    gx_loca=c_[j]#工序号
                else:#不是第一个工件
                    gx_loca=sum(self.op_num[:j])+c_[j]#工序的位置
                agv=pop3[Pi][gx_loca]#agv号
                if c[j]%2==0: #表示该工序是装载，需要判断载重,且为第一道工序
                    if a[agv]+1>AGV_weight:#超过了载重，需要更新基因的位置以及机器的位置
                        #找到上一个加工的位置的基因 然后调整到目前这个基因的位置前
                        #先更换位置，然后卸载工序
                        j=a1[agv][-1]#上一个工件号
                        for i in range(a2[agv][-1]+1,len(pop1[Pi])):
                            if pop1[Pi][i]==j:
                                b=i
                                # print("调整染色体")
                                break
                        pop1[Pi].pop(b)#删除该位置的索引的值
                        pop1[Pi].insert(Pj,j)
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


    def cross(self,pop1,pop2,pop3):
        Pop1=[]
        Pop2=[]
        Pop3=[]
        for i in range(len(pop1)//2):#遍历一半的染色体
            new_pop1 = [-1 for i in range(len(pop1[0]))]
            new_pop2 = [-1 for i in range(len(pop2[0]))]
            new_pop3 = [-1 for i in range(len(pop3[0]))]
            gj = [i  for i in range(self.fja.Jn)]
            if np.random.rand()<self.pc:
                random.shuffle(gj)#MPX交叉
                #工序交叉
                split_point=random.randint(1,len(gj)-1)
                x1=gj[:split_point]
                x2=gj[split_point:]
                for j in range(len(pop1[i])):
                    if pop1[i][j] in x1:
                        new_pop1[j]=pop1[i][j]
                site=[i for i,var in enumerate(new_pop1) if var==-1]
                e=0
                for j in range(len(pop1[i])):
                    if pop1[i+len(pop1)//2][j] in x2:
                        new_pop1[site[e]]=pop1[i+len(pop1)//2][j]
                        e=e+1
                #机器与AGV的交叉
                random_list=[random.choice([0,1]) for _ in range(sum(self.op_num))]
                for j in range(len(random_list)):
                    if random_list[j]==1:
                        new_pop2[j]=pop2[i][j]
                        new_pop3[j]=pop3[i][j]
                    else:
                        new_pop2[j] = pop2[i+len(pop1)//2][j]
                        new_pop3[j] = pop3[i+len(pop1)//2][j]
                Pop1.append(new_pop1)
                Pop2.append(new_pop2)
                Pop3.append(new_pop3)
        return Pop1,Pop2,Pop3


    def mutation(self,pop1,pop2,pop3):
        new_pop1 = copy.copy(pop1)
        new_pop2 = copy.copy(pop2)
        new_pop3 = copy.copy(pop3)
        for i in range(len(pop1)):
            if np.random.rand() <self.pm:
                random_positions=random.sample(range(len(pop1[0])),2)
                b=new_pop1[i][random_positions[0]]
                new_pop1[i][random_positions[0]]=new_pop1[i][random_positions[1]]
                new_pop1[i][random_positions[1]]=b

                random_positions1=random.sample(range(len(pop3[0])),3)
                for j in random_positions1:
                    # new_pop2[j]=random.choice[i for i in range(self.Machine) if i=new_pop2[j] and ]#不更换机器了，比较复杂
                    new_pop3[i][j]=random.choice([s for s in range(self.An) if s!=new_pop3[i][j]])
        return new_pop1,new_pop2,new_pop3


    def main(self,pop1,pop2,pop3):
        # 求适应度
        func_value1 = []
        func_value2 = []
        for i in range(len(pop1)):
            fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk,
                           epik, ea_load, ea_unload, ea_idle)
            self.Jobs, self.Machines, self.AGVs = fja.decode(pop1[i], pop2[i], pop3[i])
            v1, v2 = fja.total_value()
            func_value1.append(v1)
            func_value2.append(v2)
        non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        crowding_distance_values = []
        for j in range(0, len(non_dominated_sorted_solution)):
            crowding_distance_values.append(
                crowding_distance(func_value1[:], func_value2[:], non_dominated_sorted_solution[j][:]))  # 拥挤度计算
        # 存储层级rank与拥挤度crowd
        rank = [-1 for i in range(len(func_value1))]
        crowd = [-1 for i in range(len(func_value1))]
        for j in range(0, len(non_dominated_sorted_solution)):
            for k in range(0, len(non_dominated_sorted_solution[j])):
                rank[non_dominated_sorted_solution[j][k]] = j
                crowd[non_dominated_sorted_solution[j][k]] = crowding_distance_values[j][k]
        # 锦标赛选择个体
        new_solution = NS.Mating(len(func_value1), rank, crowd)
        new_s1 = []
        new_s2 = []
        new_s3 = []
        for n in new_solution:
            new_s1.append(pop1[n])
            new_s2.append(pop2[n])
            new_s3.append(pop3[n])
        new_pop1, new_pop2, new_pop3 = NS.cross(pop1, pop2, pop3)  # 交叉
        new_pop1 = NS.jbchrom_reset(new_pop1, new_pop3)  # 修整工序染色体
        new_pop11, new_pop21, new_pop31 = NS.mutation(pop1, pop2, pop3)  # 变异
        new_pop11 = NS.jbchrom_reset(new_pop11, new_pop31)  # 修整工序染色体
        # 交叉变异
        # 修整染色体
        char1 = new_s1 + pop1 + new_pop1 + new_pop11
        char2 = new_s2 + pop2 + new_pop2 + new_pop21
        char3 = new_s3 + pop3 + new_pop3 + new_pop31

        # 求适应度
        func_value1 = []
        func_value2 = []
        for i in range(len(char1)):
            fja = FJSP_agv(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk, epik,
                           ea_load, ea_unload, ea_idle)
            self.Jobs, self.Machines, self.AGVs = fja.decode(char1[i], char2[i], char3[i])
            v1, v2 = fja.total_value()
            func_value1.append(v1)
            func_value2.append(v2)
        non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        crowding_distance_values = []
        for j in range(0, len(non_dominated_sorted_solution)):
            crowding_distance_values.append(
                crowding_distance(func_value1[:], func_value2[:], non_dominated_sorted_solution[j][:]))  # 拥挤度计算
        # 存储层级rank与拥挤度crowd
        rank = [-1 for i in range(len(func_value1))]
        crowd = [-1 for i in range(len(func_value1))]
        for j in range(0, len(non_dominated_sorted_solution)):
            for k in range(0, len(non_dominated_sorted_solution[j])):
                rank[non_dominated_sorted_solution[j][k]] = j
                crowd[non_dominated_sorted_solution[j][k]] = crowding_distance_values[j][k]

        new_solution = []
        for i in range(0, len(non_dominated_sorted_solution)):
            non_dominated_sorted_solution2_1 = [
                index_of(non_dominated_sorted_solution[i][j], non_dominated_sorted_solution[i]) for j in
                range(0, len(non_dominated_sorted_solution[i]))]  # 返回每一层的最小值的索引，这个写的太繁琐，其实就是变成，0，1，2，3索引
            front22 = sorted(non_dominated_sorted_solution2_1, key=lambda x: crowding_distance_values[i][x])  # 根据拥挤度距离的大小进行排序，先小后打，索引需要反转
            front = [non_dominated_sorted_solution[i][front22[j]] for j in
                     range(0, len(non_dominated_sorted_solution[i]))]
            front.reverse()
            for value in front:
                new_solution.append(value)
                if (len(new_solution) == self.Pop_size):  # 其实就是根据排序和拥挤度取到足够的个体
                    break
            if (len(new_solution) == self.Pop_size):
                break
        objective_1 = [func_value1[i] for i in new_solution]
        objective_2 = [func_value2[i] for i in new_solution]
        solution_1 = [char1[i] for i in new_solution]
        solution_2 = [char2[i] for i in new_solution]
        solution_3 = [char3[i] for i in new_solution]
        return solution_1,solution_2,solution_3,objective_1,objective_2



if __name__=="__main__":
    # files = [ 11, 13, 15]
    # for file_index in range(7,16):
    #     print(f'../1_Brandimarte/BrandimarteMk{file_index}.fjs')
    #     P_time, Job_number, Machine_number, AGV_Trans, AGV_number, AGV_weight, epk, epik, ea_load, ea_unload, ea_idle = Get_Problem(
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

    NS=NSGAII(P_time, Job_number, Machine_number, AGV_Trans, AGV_number,  AGV_weight, epk, epik, ea_load,
            ea_unload, ea_idle)
    pop1, pop2, pop3 = NS.initial_population()  # 初始化
    pop1 = NS.jbchrom_reset(pop1, pop3)  # 修整工序染色体
    pareto_front_solution = []
    func1_best = []
    func2_best = []
    for g in range(GEN):
        p1, p2, p3, func_value1, func_value2 = NS.main(pop1, pop2, pop3)
        pop1 = p1
        pop2 = p2
        pop3 = p3
        non_dominated_sorted_solution = fast_non_dominated_sort(func_value1, func_value2)  # 非支配排序
        print(non_dominated_sorted_solution[0])
        non_dominated_solution = []
        for i in non_dominated_sorted_solution[0]:
            pop = []
            pop.append(pop1[i])
            pop.append(pop2[i])
            pop.append(pop3[i])
            non_dominated_solution.append(pop)
        objectives = []
        for j in non_dominated_sorted_solution[0]:
            object = []
            object.append(func_value1[j])
            object.append(func_value2[j])
            objectives.append(object)
        pareto_front_solution.append(objectives)
        ob = np.asarray(objectives)
        func1_best.append(min(ob[:, 0]))
        func2_best.append(min(ob[:, 1]))
        print(pareto_front_solution)
        print(func1_best)
    save_plot_main_data("pic_data/NSGAII/data/", pareto_front_solution,func1_best,func2_best,GEN)







