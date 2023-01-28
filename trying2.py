import random
import copy
import time
"""
INITIALIZATION PHASE
"""
import pickle
class DSMO:
    global trrr
    with open("myfile.pickle", "rb") as infile:
        path = pickle.load(infile)
    infile.close()
    global ko
    ko = time.time()
    city_dist = []
    print(path)
    f = open(path, "r")
    lines = f.readlines()
    for line in lines:

        if line == 'EOF':
            break
        newline = line[:-1].split(' ')
        for cord in newline:
            if cord != '':
                x = float(cord)
                city_dist.append(x)

    """
    Constructing the coordinate matrix
    We store the matrix as a dictonary where key is the city_id and value is a tuple containing the x and y coordinates of the city
    """
    global N
    N = len(city_dist) // 3  # N = number of cities
    global coordinates
    coordinates = {}
    for i in range(N):
        x_coordinate = city_dist[i * 3 + 1]
        y_coordinate = city_dist[i * 3 + 2]
        coordinates[city_dist[i * 3]] = (x_coordinate, y_coordinate)
    global cost
    def cost(i, j):  # calculates the cost of travelling from city i to city j
        (x_i, y_i) = coordinates[i]
        (x_j, y_j) = coordinates[j]
        c = ((x_i - x_j) ** 2 + (y_i - y_j) ** 2) ** 0.5
        return c

    global SM
    SM = []  # The list that stores all the spider monkeys
    nSM = 100  # Number of spider monkeys
    global number_of_groups
    number_of_groups = 1  # Current Number of groups
    MG = 5  # Maximum number of groups
    pr = 0.1  # Perturbation rate
    LLL = 75  # Local Leader Limit
    GLL = 75  # Global Leader Limit
    global LLLc
    LLLc = {0: 0}  # Local Leader Limit Count
    global GLLc
    GLLc = 0  # Global Leader Limit Count
    global LL
    LL = []  # List that stores the local leaders of all groups
    global GL
    GL = [0, 0]  # List that stores the group number and index of global leader

    def generate_TSP_tours(n):  # where n is the number of spider monkeys             #generates TSP tours
        temp = []
        for i in range(n):
            l = random.sample(range(1, N + 1), N)
            temp.append(l)
        SM.append(temp)

    generate_TSP_tours(nSM)
    global fitness

    def fitness(sm):
        fit = 0
        for j in range(len(sm) - 1):
            fit += cost(sm[j], sm[j + 1])
        return fit

    def SS(T1, T2):  # returns Swap sequence list with its element in a tuple e.g [(1,3),(3,4),...]
        swap_operator_list = []

        for i in range(len(T1) - 1):
            if T1[i] != T2[i]:
                first = T1[i]
                second = T2[i]

                v = T1.index(second)

                T1[i] = second
                T1[v] = first

                t = (i + 1, v + 1)

                swap_operator_list.append(t)
        return swap_operator_list

    def SS_merge(ss1, ss2):  # performs ⊗ merging operation and returns list of SS
        for i in ss2:
            if i not in ss1:
                ss1.append(i)
        return ss1

    def Cal_BasicSS(SS):  # performs Equation 7 operation and returns Basic SS list
        size = len(SS)
        # print(size)
        i = 0

        while (i != size):
            if (SS[i][1], SS[i][0]) in SS[i:]:
                SS.remove((SS[i][1], SS[i][0]))
                SS.remove(SS[i])
                size -= 2
            else:
                i += 1
        return SS

    global apply_BSS

    def apply_BSS(sm, BSS):  # Performs Equation 8 and returns updated SM

        for i in BSS:
            a = i[0] - 1
            b = i[1] - 1
            copy_sm = copy.deepcopy(sm)
            copy_sm[b], copy_sm[a] = copy_sm[a], copy_sm[b]
            if fitness(copy_sm) < fitness(sm):
                sm[b], sm[a] = sm[a], sm[b]
        return sm

    def random_subsequence(ss):
        temp = []
        for i in ss:
            u = random.randint(0, 1)
            if u == 1:
                temp.append(i)
        return temp

    iteration_counter = 1
    max_iter = 10

    LL = []
    for k in range(len(SM)):
        LLLc[k] = 0
        max_fit = 1000000000
        max_fit_index = 0
        for i in range(len(SM[k])):
            if fitness(SM[k][i]) < max_fit:
                max_fit = fitness(SM[k][i])
                max_fit_index = i
        LL.append(max_fit_index)

    GL[1] = LL[0]
    global v
    v=[]
    def dsmo_main(self):
        print(self.path)
        global LLLc
        global GLLc
        global LL
        global GL
        global number_of_groups
        global SM
        joo = time.time()
        """
        Algorithm 2.1 Starts
        """
        for k in range(len(SM)):
            for i in range(len(SM[k])):
                u = random.uniform(0, 1)
                if u >= self.pr:
                    copy_smi = copy.deepcopy(SM[k][i])

                    ss1 = self.SS(SM[k][LL[k]], SM[k][i])
                    u_ss1 = self.random_subsequence(ss1)

                    rnd = random.randint(0, len(SM[k]) - 1)
                    ss2 = self.SS(SM[k][rnd], SM[k][i])
                    u_ss2 = self.random_subsequence(ss2)

                    ssi = self.SS_merge(u_ss1, u_ss2)

                    bss_i = self.Cal_BasicSS(ssi)
                    sm_new = apply_BSS(copy_smi, bss_i)

                    if fitness(sm_new) < fitness(SM[k][i]):
                        SM[k][i] = copy.deepcopy(sm_new)

        for k in range(len(SM)):
            for i in range(len(SM[k])):
                copy_smi = copy.deepcopy(SM[k][i])

                ss1 = self.SS(SM[GL[0]][GL[1]], SM[k][i])
                u_ss1 = self.random_subsequence(ss1)

                rnd = random.randint(0, len(SM[k]) - 1)
                ss2 = self.SS(SM[k][rnd], SM[k][i])
                u_ss2 = self.random_subsequence(ss2)

                ssi = self.SS_merge(u_ss1, u_ss2)

                bss_i = self.Cal_BasicSS(ssi)
                sm_new = apply_BSS(copy_smi, bss_i)

                if fitness(sm_new) < fitness(SM[k][i]):
                    SM[k][i] = copy.deepcopy(sm_new)

        """
        Algorithm 2.2 Starts
        """
        best_LL_group = 0
        best_LL_index = 0
        best_LL_fitness = fitness(SM[GL[0]][GL[1]])
        for k in range(len(SM)):
            LLk_fitness = fitness(SM[k][LL[k]])
            new_LL_maximum = LLk_fitness
            new_LL_index = 0
            for i in range(len(SM[k])):
                if i == LL[k]:
                    continue
                fit = fitness(SM[k][i])
                if fit < new_LL_maximum:
                    new_LL_maximum = fit
                    new_LL_index = i

            if new_LL_index == 0:
                LLLc[k] += 1
            else:
                LLLc[k] = 0
                LL[k] = new_LL_index
                self.LLk_fitness = new_LL_maximum

            if LLk_fitness < best_LL_fitness:
                best_LL_group = k
                best_LL_index = LL[k]
                best_LL_fitness = LLk_fitness

        if best_LL_group != GL[0] or best_LL_index != GL[1]:
            GL[0] = best_LL_group
            GL[1] = best_LL_index
            self.GLLc = 0

        else:
            self.GLLc += 1

        """
        Algorithm 2.3 Starts
        """
        for k in range(len(SM)):
            if LLLc[k] > self.LLL:
                LLLc[k] = 0
                for i in range(len(SM[k])):
                    u = random.uniform(0, 1)
                    if u >= self.pr:
                        l = random.sample(range(1, N + 1), N)
                        SM[k][i] = l
                    else:
                        copy_smi = copy.deepcopy(SM[k][i])

                        ss1 = self.SS(SM[k][LL[k]], SM[k][i])
                        u_ss1 = self.random_subsequence(ss1)

                        ss2 = self.SS(SM[GL[0]][GL[1]], SM[k][i])
                        u_ss2 = self.random_subsequence(ss1)

                        sm_new = apply_BSS(copy_smi, u_ss2)
                        sm_new2 = apply_BSS(sm_new, u_ss1)
                        SM[k][i] = sm_new2

        if self.GLLc > self.GLL:
            GLLc = 0

            if self.number_of_groups <self. MG:
                lis = []
                for k in range(len(SM)):
                    for i in range(len(SM[k])):
                        lis.append(SM[k][i])
                self.number_of_groups += 1
                SM = []
                k = 0
                x = self.nSM // self.number_of_groups
                temp = []
                for i in range(self.nSM):
                    if (i % x == 0) and (k != self.number_of_groups - 1) and (i != 0):
                        SM.append(temp)
                        temp = []
                        k += 1
                    temp.append(lis[i])
                SM.append(temp)

            else:
                lis = []
                for k in range(len(SM)):
                    for i in range(len(SM[k])):
                        lis.append(SM[k][i])
                SM = []
                SM.append(lis)
                number_of_groups = 1

            LL = []
            LLLc = {}
            for k in range(len(SM)):
                LLLc[k] = 0
                max_fit = 10000000000000
                max_fit_index = 0
                for i in range(len(SM[k])):
                    if fitness(SM[k][i]) < max_fit:
                        max_fit = fitness(SM[k][i])
                        max_fit_index = i
                LL.append(max_fit_index)

            GL = [0, 0]
            for k in range(len(SM)):
                if fitness(SM[k][LL[k]]) < fitness(SM[GL[0]][GL[1]]):
                    GL[0] = k
                    GL[1] = LL[k]
        bet=fitness(SM[GL[0]][GL[1]])
        print("iter num {} ".format(self.iteration_counter),bet , "\n", (time.time() - joo))


        v.append(bet)
        print(min(v), "\ntotal time: ", (time.time() - ko), " s")


        if self.iteration_counter==11 or self.iteration_counter==101:
            self.trrr=min(v)
        else:
            self.trrr = bet
        self.iteration_counter += 1
    def sed(self):
        return self.trrr
    """f.write(str(run))
    f.write(" ")
    f.write(str(fitness(SM[GL[0]][GL[1]])))
    f.write("\n")
    f.write(str(SM[GL[0]][GL[1]]))
    f.write("\n")"""
    # f.close()








