def homos(S,f):
    """String homomorphism wrt lambda f
    homos("abcd",hm) --> 'bcde' where hm = lambda x: chr( (ord(x)+1) % 256 )
    Returns a homomorphism of map(S,f)
    """
    return "".join(map(f,S))

def dotsan_map(x):
        """Remove brackets and commas and separate elements with '_'
        Does not return a homomorphism
        """
        if x in { "{", " ", "'", "}" }:
            return ""
        elif x == ",":
            return "_"
        else:
            return x

def dot_san_str(S):
        """Make dot like strings which are in set of states notation."""
        return homos(S, dotsan_map)

def prDotHeader(fl):
        """Prints the header of the .DOT file for DFAs"""
        fl.write(r'digraph G {' + "\n")
        fl.write(r'/* Defaults */' + "\n")
        fl.write(r' fontsize = 12;' + "\n")
        fl.write(r' ratio = compress; ' + "\n")
        fl.write(r' rankdir=LR; ' + "\n")
        fl.write(r'/* Bounding box */' + "\n")
        fl.write(r' size = "4,4";' + "\n")

def prNonFinalNodeName(fl, q):
        """Prints non accepting nodes"""
        fl.write(dot_san_str(q) + r'[shape=circle, peripheries=1];' + "\n")

def prFinalNodeName(fl, q):
        """Prints accepting nodes """
        # Could write like print (q, r'[shape=circle, peripheries=2];', file=fl)
        # But am documenting use of trailing comma to suppress \n . In Python3 we supply end = ''
        fl.write(dot_san_str(q)) # end with no CR
        fl.write(r'[shape=circle, peripheries=2];' + "\n") # end with a CR

def prOrientation(fl):
        """Prints the horizontal/vertical orientation of the file
        """
        fl.write(r'/* Orientation */' + "\n")
        fl.write(r'orientation = landscape;' + "\n")

def prEdges_w_bh(fl, D):
        """Prints edges;
        Edges that lead to a black hole are included
        """
        fl.write(r'/* The graph itself */' + "\n")
        fl.write(r'"" -> ', dot_san_str(D["q0"]), ";" + "\n")
    
        for QcQ in D["Delta"].items():
            fl.write(dot_san_str(QcQ[0][0])+ r' -> '+
                dot_san_str(QcQ[1]) + r'[label="' + dot_san_str(QcQ[0][1]) + r'"];' + "\n")

def prEdges(fl, D):
        """Prints edges;
        Edges that lead to a black hole are not included
        """
        fl.write(r'/* The graph itself */' + "\n")
        fl.write(r'"" -> ' + dot_san_str(D["q0"]) + ";" + "\n")
        for QcQ in D["Delta"].items():
            if (((QcQ[0][0]) != "BH") & (QcQ[1] != "BH")):
                fl.write(dot_san_str(QcQ[0][0]) + r' -> ' +
                    dot_san_str(QcQ[1]) + r'[label="' + dot_san_str(QcQ[0][1]) + r'"];' + "\n")

def prClosing(fl):
        """Prints the closing text of the .DOT file
        """
        fl.write(r'/* Unix command: dot -Tps exdfa.dot >! exdfa.ps */' + "\n")
        fl.write(r"/* For further details, see the `dot' manual */" + "\n")
        fl.write(r"}" + "\n")

def prNodeDefs_w_bh(fl, D):
        """Prints all of the nodes along with where their arrow leads
        Single Circles represent non-accepting nodes
        Black Holes are included
        """

        fl.write(r'/* Node definitions */' + "\n")
        fl.write(r' "" [shape=plaintext];' + "\n") # Start state arrow is from "" to I
        # All non-accepts are single circles
        for q in D["Q"] - D["F"]:
            prNonFinalNodeName(fl, q)
        for q in D["F"]:
            prFinalNodeName(fl, q)

def prNodeDefs(fl, D):
        """Prints all of the nodes along with where their arrow leads
        Single Circles represent non-accepting nodes
        Black Holes are not included
        """
        fl.write(r'/* Node definitions */' + "\n")
        fl.write(r' "" [shape=plaintext];' + "\n") # Start state arrow is from "" to I
        # All non-accepts are single circles
        for q in D["Q"] - D["F"]:
            if (q != "BH"):
                prNonFinalNodeName(fl, q)
        for q in D["F"]:
            prFinalNodeName(fl, q)

def dot_dfa(D, fname):
        """Generate a dot file with the automaton in it. Run the dot file through
        dot and generate a ps file."""
        fl = open(fname, 'w')
        #-- digraph decl
        prDotHeader(fl)
        #-- node names and how to draw them
        prNodeDefs(fl, D)
        #-- orientation - now landscape
        prOrientation(fl)
        #-- edges
        prEdges(fl, D)
        #-- closing
        prClosing(fl)

if __name__ == "__main__":
    Q1 = {'S0','S1'}
    Sigma1 = {'a','b'}
    Delta1 = {('S0', 'a'): 'S0', ('S1', 'a'): 'S0', ('S1', 'b'): 'S1', ('S0', 'b'): 'S1'}
    q01 = 'S0'
    fl = {'S1'}
    dot_dfa(mk_dfa(Q1,Sigma1,Delta1,q01,fl),'DotDFA.dot')

    help(homos)
    help(dotsan_map)
    help(dot_san_str)
    help(prDotHeader)
    help(prNonFinalNodeName)
    help(prFinalNodeName)
    help(prOrientation)
    help(prEdges_w_bh)
    help(prEdges)
    help(prClosing)
    help(prNodeDefs_w_bh)
    help(prNodeDefs)
    help(dot_dfa)