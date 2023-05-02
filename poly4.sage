print("ver 2022.10.14 (poly2)")

#several functions now form part of Poly and are therefore not included here


from functools import cmp_to_key

def face_comp(x,y):
    nx=len([e.rev.face for e in x._edges if not e.is_border])
    ny=len([e.rev.face for e in y._edges if not e.is_border])
    if nx>ny:
        return int(-1)
    if nx==ny and x.label<y.label:
        return int(-1)
    if nx==ny and x.label==y.label:
        return int(0)
    return int(1)

def transversal_region_bad(P,H):
    #P is expected to be a Poly and H a subgroup of P.group()
    #The goal is to get a connected fundamental domain for X/H consisting of a subset of the faces of the fundamental domain for X.
    #This corresponds to finding right coset representatives of H
    gens=P.generators
    G=P.group()
    S=P.faces
    T=[S[0]]
    S=[F for F in S if not F.g in H] #faces not yet represented in T
    TNeighbors=[e.rev.face for k in T for e in k._edges if e.rev.face in S and not e.is_border]  #faces connected to T
    while len(S)>0:
        f=TNeighbors[0]
        T.append(f)
        for e in f._edges:
            f_n=e.rev.face
            if not e.is_border and f_n in S and not f_n in T and not f_n in TNeighbors:
                TNeighbors.append(f_n)
        S=[F for F in S if not F.g*(f.g).inverse() in H]
        TNeighbors=[F for F in TNeighbors if F in S]
    if len(T)*H.order()==G.order():
        #return GeodesicRegion(deepcopy(T),gens)
        picture = P.select_faces_by_words([k.word for k in T])
        return GeodesicRegion(picture,gens)
    return "ERROR"

def transversal_region(P,H):
    #P is expected to be a Poly and H a subgroup of P.group()
    gens=P.generators
    G=P.group()
    attempt=0
    contador=0
    S=set(G).difference(set(H))
    T=[P.faces[0]]
    TNeighbors=[e.rev.face for k in T for e in k._edges if not e.is_border]
    TNeighbors.sort(key=cmp_to_key(face_comp))
    progress=True
    while len(S)>0 and progress:
        progress=False
        contador+=1
        while len(TNeighbors)>0:
            f=TNeighbors[0]
            x=f.g
            if x in S:
                progress=True
                T.append(f)
                for e in f._edges:
                    if not e.is_border and not e.rev.face in T and not e.rev.face in TNeighbors:
                        TNeighbors.append(e.rev.face)
                TNeighbors.remove(f)
                TNeighbors.sort(key=cmp_to_key(face_comp))
                for h in H:
                    S.discard(G(h)*x)
                break
            else:
                TNeighbors.remove(f)
    if len(T)*H.order()==G.order():
        return GeodesicRegion(T,gens)
    return GeodesicRegion(T,gens),S

def induced_generating_vector(G,v,H):
    A=[k.order() for k in v]  #signature
    gX=1+(0-1)*G.order()+G.order()/2*sum(1-1/k for k in A) #genus of curve on which G acts
    P=Poly(G,v)
    L1=transversal_region(P,H) #corresponding to right cosets of H in G
    if L1=="ERROR":
        print(G.id(),H.id(),"error no entendido")
        return -1,-1,-1
    L=[[k.g.inverse()*j for j in H]for k in L1]  #left cosets of H in G
    GH=dict()#lookup table to find which coset in L contains a given g in G
    w=[]
    wx=[]
    for c in range(len(L)):
        for g in L[c]:
            GH[g]=c  #assigns a number corresponding to the coset in which g is
    B=L1.border()
    expg=L1._expanded_generators
    available=[[True for j in expg] for k in L]
    for f in L1:  #this is a fix related to the split elliptic point
        ww=f._edges[0].face.word
        f._edges[0].a().w=ww*P.fundamental_polygon.letters()[-1]*ww.inverse()
        f._edges[0].a().label=f._edges[0].face.g
    for e in B:
        if e.a().w!=None:  #selecting edges starting at elliptic points
            eg=(e.a().w)(v)  #image in G of fuchsian group element fixing elliptic point
            f=GH[e.face.g.inverse()]  #position in L of coset
            edir=e.direction if e.direction%2 else -2  #split elliptic point has to be treated separately
            if edir==-2:
                f=GH[e.a().label.inverse()]
            a=expg[edir] 
            z=GH[a*L[f][0]]
            k=1
            while z!=f:
                available[z][edir]=False
                k+=1
                z=GH[a*L[z][0]]
            if available[f][edir] and eg.order()>k:
                available[f][edir]=False
                w.append(eg^k)
                wx.append((e.a().w)^k)
    A=[k.order() for k in w]
    gH=(gX-1)/H.order()-1/2*sum(1-1/k for k in A)+1  #genus of X/H
    return w,wx,gH,A

