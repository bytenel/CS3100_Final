from dfa import *
from praut import *

        
def minDFA(D, state_name_mode):
    """Given a DFA D, go through the state minimization algorithm.
    state_name_mode is 'verbose' or 'succinct', producing two variants, as you can guess.
    """
    #
    # Implement the code for minDFA
    #

    

#=================================================================

def mk_DFAFig1011():
    Q = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}
    Sigma = {'a','b'}
    Delta = {('S1', 'a'): 'S2',
             ('S1', 'b'): 'S3',
             ('S2', 'a'): 'S4',
             ('S3', 'a'): 'S5',
             ('S2', 'b'): 'S5',
             ('S3', 'b'): 'S4',
             ('S5', 'a'): 'S6',
             ('S5', 'b'): 'S6',
             ('S4', 'a'): 'S6',
             ('S4', 'b'): 'S6',
             ('S6', 'a'): 'S6',
             ('S6', 'b'): 'S6'}
    q0 = 'S1'
    F = {'S2', 'S3', 'S6'}
    return mk_dfa(Q, Sigma, Delta, q0, F)


def mk_tree2DFA():
    Q = {'S','Sa', 'Sb', 'Sc', 'Saa', 'Scc'}
    Sigma = {'a','b','c'}
    Delta = {('S', 'a'): 'Sa',
             ('S', 'b'): 'Sb',
             ('S', 'c'): 'Sc',
             ('Sa', 'a'): 'Saa',
             ('Sc', 'c'): 'Scc'}
    q0 = 'S'
    F = {'Saa', 'Scc'}
    return mktot(mkp_dfa(Q, Sigma, Delta, q0, F))

def mk_Qtt():
    Q = {'S','Sa', 'Sb', 'Sc', 'Saa', 'Sab', 'Sac'}
    Sigma = {'a','b','c'}
    Delta = {('S', 'a'): 'Sa',
             ('S', 'b'): 'Sb',
             ('S', 'c'): 'Sc',
             ('Sa', 'a'): 'Saa',
             ('Sa', 'b'): 'Sab',
             ('Sa', 'c'): 'Sac' }
    q0 = 'S'
    F = {'Saa', 'Sab', 'Sac'}
    return mktot(mkp_dfa(Q, Sigma, Delta, q0, F))


def mk_tree3DFA():
    Q = {'S','Sa', 'Sb', 'Sc', 'Saa', 'Sab', 'Sac', 'Sba', 'Sbb', 'Sbc', 'Sca', 'Scb', 'Scc',
          'Saaa', 'Sccc'}
    Sigma = {'a','b','c'}
    Delta = {('S', 'a'): 'Sa',
             ('S', 'b'): 'Sb',
             ('S', 'c'): 'Sc',
             ('Sa', 'a'): 'Saa',
             ('Sa', 'b'): 'Sab',
             ('Sa', 'c'): 'Sac',
             ('Sb', 'a'): 'Sba',
             ('Sb', 'b'): 'Sbb',
             ('Sb', 'c'): 'Sbc',
             ('Sc', 'a'): 'Sca',
             ('Sc', 'b'): 'Scb',
             ('Sc', 'c'): 'Scc',
             ('Saa', 'a'): 'Saaa',
             ('Scc', 'c'): 'Sccc' }
    q0 = 'S'
    F = {'Saaa', 'Sccc'}
    return mktot(mkp_dfa(Q, Sigma, Delta, q0, F))

minDTree3 = minDFA(mk_tree3DFA(), 'verbose')


#---min DFA
# ==> Enable if you wish : dot_dfa(minDTree3, "minDTree3.dot")

#=================================================================






    
    

    
