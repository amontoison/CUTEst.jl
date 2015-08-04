#!/usr/bin/env python

import re
import sys

foos = [ "ccfg", "ccfsg", "ccifg", "ccifsg", "cdh", "cfn", "cgr", "cgrdh",
        "chcprod", "chprod", "cidh", "cish", "cjprod", "clfg", "cofg", "cofsg",
        "csgr", "csgrsh", "csh", "cshc" ]

trip_comp = [ "j_var", "j_fun", "h_row", "h_col" ]
triplet = { "Jx": "J", "Wx": "W" }
trip_resp = { "Jx": "j_var", "Wx": "h_row" }

ignore = trip_comp + [ "nnzj", "nnzh", "nnzg", "nnzgci", "gci_var", "g_var" ]

translate = {
    "n": "nlp.meta.nvar",
    "m": "nlp.meta.ncon",
    "lgci": "nlp.meta.nvar",
    "x": "x0",
    "c": "cx",
    "f": "fx",
    "g": "gx",
    "grad": "true",
    "goth": "false",
    "gotj": "false",
    "cjac": "Jx",
    "icon": "j",
    "iprob": "j",
    "j_val": "Jx",
    "h_val": "Wx",
    "jtrans": "false",
    "lcjac1": "nlp.meta.ncon",
    "lcjac2": "nlp.meta.nvar",
    "lg": "nlp.meta.nvar",
    "y": "y0",
    "grlagf": "false",
    "lj": "nlp.meta.nnzj+nlp.meta.nvar",
    "lj1": "nlp.meta.ncon",
    "lj2": "nlp.meta.nvar",
    "lh1": "nlp.meta.nvar",
    "lh": "nlp.meta.nnzh",
    "lvector": "nlp.meta.nvar",
    "lresult": "nlp.meta.ncon",
    "vector": "ones(nlp.meta.nvar)" }

hs = {
    "fx": "f(x0)",
    "cx": "c(x0)",
    "ci": "c(x0)[j]",
    "h": "(W(x0,[i == j ? 1.0 : 0.0 for i = 1:nlp.meta.ncon])-H(x0))",
    "gx": "g(x0)",
    "Jx": "J(x0)",
    "gci": "J(x0)[j,:]",
    "g_val": "g(x0)[g_var]",
    "gci_val": "J(x0)[j,gci_var]",
    "Wx": "W(x0,y0)" }

special = {
    "cshc": { "W": "W(x0,y0)-H(x0)"},
    "cish": { "W": "W(x0,[i == j ? 1.0 : 0.0 for i = 1:nlp.meta.ncon])-H(x0)" },
    "clfg": { "fx": "f(x0)+dot(y0,c(x0))", "gx": "g(x0)+J(x0)'*y0" }
    }

results = {
    "chprod": "W(x0,y0)*v",
    "chcprod": "(W(x0,y0)-H(x0))*v",
    "cjprod": "J(x0)*v" }

sizeof = {
    "gx": "nlp.meta.nvar",
    "cx": "nlp.meta.ncon",
    "Jx": "nlp.meta.ncon, nlp.meta.nvar",
    "Wx": "nlp.meta.nvar, nlp.meta.nvar",
    "j_var": "Int, nlp.meta.nnzj+nlp.meta.nvar",
    "j_fun": "Int, nlp.meta.nnzj+nlp.meta.nvar",
    "h_row": "Int, nlp.meta.nnzh",
    "h_col": "Int, nlp.meta.nnzh",
    "gci_val": "nlp.meta.nvar",
    "gci_var": "Int, nlp.meta.nvar",
    "g_val": "nlp.meta.nvar",
    "g_var": "Int, nlp.meta.nvar",
    "h": "nlp.meta.nvar, nlp.meta.nvar" }


sizeofsp = {
    "Jx": "nlp.meta.nnzj+nlp.meta.nvar",
    "Wx": "nlp.meta.nnzh" }

multiples = {
    "icon": "nlp.meta.ncon",
    "iprob": "nlp.meta.ncon" }

zero_index = [ "csgr", "csgrsh" ]

def addTriplet(trip, spc, fname):
    str = spc+"{}_val = copy({}x)\n".format(trip.lower(), trip)
    if trip == "J":
        str += spc+"Jx = zeros(nlp.meta.ncon, nlp.meta.nvar)\n"
        i = "j_fun"
        j = "j_var"
        z = "j"
    else:
        str += spc+"{}x = zeros(nlp.meta.nvar, nlp.meta.nvar)\n".format(trip)
        i = "h_row"
        j = "h_col"
        z = "h"
    str += spc+"for k = 1:nnz{}\n".format(z)
    if fname in zero_index:
        str += spc+"  {}[k] == 0 && continue\n".format(i)
    str += spc+"  {}x[{}[k],{}[k]] = {}_val[k]\n".format(trip, i, j, trip.lower())
    if trip != "J":
        str += spc+"  {}x[{}[k],{}[k]] = {}_val[k]\n".format(trip, j, i, trip.lower())
    str += spc+"end\n"
    if fname in special and trip in special[fname]:
        str += spc+"@test_approx_eq_eps {}x {} 1e-8\n".format(trip, \
                special[fname][trip])
    else:
        str += spc+"@test_approx_eq_eps {}x {} 1e-8\n".format(trip, \
                hs[trip+"x"])
    return str

def generate_test_for_function (foo):
    match = ""
    for mult in multiples.keys():
        if mult in foo:
            match = mult
            break
    if match != "":
        head = "for j = 1:" + multiples[match] + "\n"
        spc = "  "
    else:
        spc = ""
        head = ''
    str = ""
    fname = foo.split("(")[0].strip()
    inplace = fname[-1] == "!"
    cint = "Cint" in foo.split(")")[0]
    inputs = []
    for x in foo[foo.find("(")+1:foo.find(")")].split(","):

        if "libname" in x:
            inputs.append("nlp.libname")
        else:
            x = x[0:x.find(":")].strip()
            if x not in ["1", "2", "3"]:
                x = translate[x] if x in translate else x
                inputs.append(x)
    outputs = []
    for x in foo[foo.find("return")+6:-6].split(","):
        if "[" in x:
            x = x[0:x.find("[")]
        x = x.strip()
        if len(x) == 0:
            continue
        x = translate[x] if x in translate else x
        outputs.append(x)

    str += spc
    if len(outputs) > 0:
        str += ', '.join(outputs) + " = "
    str += "{}({})\n".format(fname, ', '.join(inputs))

    if inplace:
        fname = fname[0:len(fname)-1]

    for x in outputs:
        if x in triplet and trip_resp[x] in outputs:
            str += addTriplet(triplet[x], spc, fname)
        elif x not in ignore:
            if x == "result":
                str += spc+"@test_approx_eq_eps {} {} 1e-8\n".format(x, results[fname])
            elif fname in special and x in special[fname]:
                str += spc+"@test_approx_eq_eps {} {} 1e-8\n".format(x, \
                        special[fname][x])
            else:
                str += spc+"@test_approx_eq_eps {} {} 1e-8\n".format(x, hs[x])
    if inplace:
        for x in inputs:
            if x == "result":
                str = spc+"{} = zeros({})\n".format(x,results[fname]) + str
            if x in sizeof:
                if x in triplet and trip_resp[x] in inputs:
                    str = spc+"{} = zeros({})\n".format(x,sizeofsp[x]) + str
                    str += addTriplet(triplet[x], spc, fname)
                elif fname in special and x in special[fname]:
                    str += spc+"@test_approx_eq_eps {} {} 1e-8\n".format(x, \
                            special[fname][x])
                else:
                    arg = sizeof[x]
                    if cint and "Int" in sizeof[x]:
                        arg = arg.replace("Int", "Cint")
                    str = spc+"{} = zeros({})\n".format(x,arg) + str
                    if x not in ignore:
                        str += spc+"@test_approx_eq_eps {} {} 1e-8\n".format(x, hs[x])
    if match != "":
        str += "end\n"

    return head+str+"\n"

filename = "src/specialized_interface.jl"
content = ''.join(open(filename, "r").readlines()).split("function")[1:]

selection = []
for x in content:
    if any([foo+"(" in x for foo in foos]):
        selection.append(x)

with open("test/test_specialized.jl","w") as f:
    f.write('println("\\nTesting the Specialized interface\\n")\n\n')
    f.write('v = ones(nlp.meta.nvar)\n')
    for x in selection:
        f.write(generate_test_for_function(x))