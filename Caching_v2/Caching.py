import itertools
import math    
import collections

def findsubsets(S,m):
    return list(itertools.combinations(S, m))

def nCr(n,r):
    f = math.factorial
    return int(f(n) / f(r) / f(n-r))

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

if __name__ == "__main__": 
    N = 3 # No. of files
    file = ['A','B','C','D','E']
    K = 6 # No. of users
    t = 3
    
    print('========================= Transmission ==============================')
    community = findsubsets(range(1, K+1), t+1) # all subsets of size 't+1'
    print(community)
    num_users_in_com = t+1
    demand = ['A', 'A', 'A', 'A', 'B', 'C']
    # delete any community that is redundant:
    leader = []
    for ii in file: # find the leader 
        if ii in demand:
           leader.append(demand.index(ii))
    leader = [x + 1 for x in leader] # leader index (start from 1 to K)
    print('The leaders are:', leader)
    card_of_leader = len(leader) # No. of leaders
    card_of_B = t + 1 + card_of_leader # cardinality of set B
    users = list(range(1, K+1))
    for ii in leader:
        users.remove(ii)
    if card_of_B >= len(users):
        Y_B_exclude_V = findsubsets(users, card_of_B - card_of_leader)
    print('The delated Y are:', 'Y_', ''.join(str(v) for v in Y_B_exclude_V))
    for ii in Y_B_exclude_V: # remove the redandant linear comb
        community.remove(ii)
    num_community = len(community) # No. of communities left
    print('The original No. of linear combinations is:',nCr(K,t+1))
    print('The remaining No. of linear combinations is:', num_community)
    if K-card_of_leader >= t+1:
       print('Matches with the theoretical result?:', num_community == (nCr(K,t+1)-nCr(K-card_of_leader,t+1)))
    else:
       print('Matches with the theoretical result?:', num_community == (nCr(K,t+1)))       
    num_of_x = [] # No. of variables corresponds to each linear comb
    linearCombBreak = [] # all the linear combs
    for iCom in range(num_community): # each community 
        linearComb = ''
        for iUser in range(num_users_in_com): # each user in a community
            users = list(community[iCom])
            user = users[iUser]
            del users[iUser]
            if iUser < num_users_in_com - 1:
                linearComb = linearComb + demand[user-1] + '_{' + ''.join(str(v) for v in users) + '} + ' # find this user's demand 
            else:
                linearComb = linearComb + demand[user-1] + '_{' + ''.join(str(v) for v in users) + '} \\\\'# find this user's demand 
        print(linearComb)
        # print('Break into:')
        x = 0
        for ii in file: # break up each linear combination into combi.  within the same file
            pos = find(linearComb, ii)
            if len(pos) != 0:
                linearCombBreak.append(linearComb[pos[0]:pos[-1]+1+2+t+1])
                # print(linearComb[pos[0]:pos[-1]+1+2+t+1])
                x += 1
        num_of_x.append(x)
        print('\n')
    # print(linearCombBreak)
    print('Num of variables x to create:', num_of_x)
    duplicates = [item for item, count in collections.Counter(linearCombBreak).items() if count > 1]
    if len(duplicates) != 0:
        print('There are duplicates, introduce extra constraints:')
        for ii in duplicates:
            pos = find(linearCombBreak, ii) # position that duplicates exist
            pos = [x+1 for x in pos] # indice starts from 1, instead of 0
            line = ''
            for jj in pos:
                line = line + 'x' + str(jj) + ',' 
            print('There is a constraint between', line)
    else:
         print('There are not duplicates.')
    # print([item for item, count in collections.Counter(linearCombBreak).items() if count > 1]) # find the duplicates of linear comb after break up
    

    '''
    print('========================= Cache ==============================')
    t_subsets = findsubsets(range(1,K+1), t)
    for iUser in range(K): # derive what's in each user's cache
        cache = list()
        iUser = iUser + 1 #index to user index
        print('Symbols in user', iUser, ':')
        for iSubset in range(len(t_subsets)):
            subset = list(t_subsets[iSubset])
            if iUser in subset: # if this subset includes this user
                for ifile in range(len(file)):
                    cache.append(file[ifile])
                    cache.append('_{')
                    cache.append(''.join(str(v) for v in subset))
                    cache.append('},')
        print((''.join(str(n) for n in cache))[:-1] + ' \\\\ \n')        
        '''