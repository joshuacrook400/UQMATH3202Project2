from gurobipy import *
import pandas
import math

'''
Joshua Crooks (s46974408) solution to MATH3202 Linear Programming Assignment
Date:  20/03/2023
'''

###Data###

nodes = pandas.read_csv("A1nodes2.csv")
grid = pandas.read_csv("A1pipelines.csv")

Months = ['D0','D1', 'D2', 'D3', 'D4', 'D5', 
          'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'D13']

T = range(len(Months))
N = range(len(nodes['Node']))

#matrix of demands
demand = []
for t in T:
    # List of demand for nodes for one day
    x = [nodes[Months[t]][j] for j in N]
    demand.append(x) 

# capacity and cost for each supplier node
capacity = { 0: 340, 4: 773, 15: 761, 37:553 }
cost = { 0: 67, 4: 90, 15: 84, 37:69 }

# E[i,j] gives the distance (km) from node i to node j
E = {}
for j in range(len(grid['Pipeline'])):
    n1 = grid['Node1'][j]
    n2 = grid['Node2'][j]
    distance = math.hypot(nodes['X'][n1]-nodes['X'][n2],nodes['Y'][n1]-nodes['Y'][n2])
    E[n1,n2] = distance

#creating gurobi model
m = Model("PacificParadiseGas")   

###Variables###

#gas produced at each node each day
X = {(n,t): m.addVar() for n in N for t in T}

#gas flowing in each pipe each day
Y = {(e,t): m.addVar() for e in E for t in T}

#imbalances of each pipe. Adds an extra day to allow for a deficit on the first day
B = {(e,t): m.addVar(lb = -10^8) for e in E for t in T}

#Variable used to define absolute value of B used for calculating costs on negative imbalances.
Z = {(e,t): m.addVar(lb = -10^8) for e in E for t in T}


###Objective###

#cost of buying gas + cost of transporting gas +cost of imbalances
m.setObjective(quicksum(cost[n]*X[n,t] for n in N for t in T if n in cost) 
               + quicksum(0.01 * Y[e,t] * E[e] for e in E for t in T)
               + quicksum(0.1*Z[e,t] for e in E for t in T)
               , GRB.MINIMIZE)    


###Constraints###
#Constraint logging
dailyFlowConstr = {}
dailyCapConstr = {}
imbConstr = {}

#daily constraints 
for t in T:
    
    #creating absolute value of the imbalance
    for e in E:
        m.addConstr(Z[e,t] >= B[e,t])
        m.addConstr(Z[e,t] >= -B[e,t])

    # Supplier capacity
    for n in N:
        if n in capacity:
            dailyCapConstr[(n,t)] = m.addConstr(X[n,t] <= capacity[n])
        else:
            m.addConstr(X[n,t] <= 0)

    #pipeline capacity
    for e in E:
        dailyFlowConstr[(e,t)] = m.addConstr(Y[e,t] <= 426)
        
        
        #imbalance constraint, only used for sensitivity analysis
        imbConstr[(e,t)] = m.addConstr(B[e,t] <= 426 )
        #m.addConstr(B[e,t] >= -426 )
        
    # Flow Balance constraint
    for n in N:
        #final timestep
        if t == T[-1]:
            m.addConstr( X[n,t] + quicksum(Y[e,t] + B[e,t] for e in E if e[1] == n) == 
            quicksum(Y[e,t] for e in E if e[0] == n) + demand[t][n])
        else:
            m.addConstr( X[n,t] + quicksum(Y[e,t] + B[e,t] for e in E if e[1] == n) == 
            quicksum(Y[e,t] + B[e,t+1] for e in E if e[0] == n) + demand[t][n])

suplimitConstr = {}
#fortnite supplier limit
for n in N:
    if n in cost:
        suplimitConstr[n] = m.addConstr(quicksum(X[n,t] for t in T) <= 7202) 

#Pipe imbalance net 0 
for e in E:
    m.addConstr(quicksum(B[e,t] for t in T) == 0)

#Solve the problem    
m.optimize()

print("Total cost=",round(m.objval,2))
print()
for t in T:
    print('Day',t+1,'gas',[round(X[n,t].x) for n in N if n in capacity])

print('###Sensitivity Analysis###')

print('#Suppliers#')
for n in N:
    if n in capacity:
        print('Node:',n
              ,' Total Gas:',sum(X[n,t].x for t in T)
              ,'Total Reduced Cost:',sum(X[n,t].RC for t in T)
              ,'Supplier 2 week limit dual value:', round(suplimitConstr[n].Pi,4)
              ,'Total daily supplier cap dual value ', sum(dailyCapConstr[(n,t)].Pi for t in T))
        
print('#Pipes#')
i = 1
for e in E:
    if sum(dailyFlowConstr[(e,t)].Pi for t in T) != 0:
        print('Pipe:',i
            ,'Total gas flow through:',round(sum(Y[e,t].x for t in T))
            , 'Total Reduced Cost:', round(sum(Y[e,t].RC for t in T),3)
            ,'Total Pipe flow cap dual value:',round(sum(dailyFlowConstr[(e,t)].Pi for t in T),4)
            ,'Slack:',round(sum(dailyFlowConstr[(e,t)].Slack for t in T),4)
            ,'bounds low:', round(sum(dailyFlowConstr[(e,t)].SARHSLow for t in T)))
    i = i + 1

print('#Large Imbalances#')
for t in T:
    i = 1
    for e in E:
        if abs(B[e,t].x) > 25:
            print('Day:',t+1
                  ,'Pipe:',i
                  ,'Nodes:',e
                  ,'Imbalance:',B[e,t].x
                  ,'SAUB bounds:',[B[e,t].SAUBLow, B[e,t].SAUBUp]
                  ,'Dual value:',imbConstr[(e,t)].Pi
                  ,'SARHSLow:',imbConstr[(e,t)].SARHSLow
                )
        i = i +1
 

