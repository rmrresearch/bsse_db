import pickle
from pathlib import Path
from open_pickle_data import open_pickle

def twobody_bsse(pickle_file,molecule):
    

    raw_data = open_pickle(pickle_file)
    two_body_no_vmfc = raw_data['2b_no_vmfc']
    two_body_vmfc = raw_data['2b_vmfc']
    #--
    two_body_bsse={}
    two_body_bsse[molecule] = {}
    for cluster in two_body_no_vmfc[molecule]:
        two_body_bsse[molecule][cluster] = {}
        #--
        for motion in two_body_no_vmfc[molecule][cluster]:
            two_body_bsse[molecule][cluster][motion] = {}
            #--
            for config in two_body_no_vmfc[molecule][cluster][motion]:
                two_body_bsse[molecule][cluster][motion][config]={}
                #--
                for method in two_body_no_vmfc[molecule][cluster][motion][config]:
                    two_body_bsse[molecule][cluster][motion][config][method] = {}
                    #--
                    for basis in two_body_no_vmfc[molecule][cluster][motion][config][method]:
                        two_body_bsse[molecule][cluster][motion][config][method][basis] = {}
                        #--
                        DE_no_vmfc = two_body_no_vmfc[molecule][cluster][motion][config][method][basis]['energy']
                        DE_vmfc =  two_body_vmfc[molecule][cluster][motion][config][method][basis]['energy']
                        #--
                        two_body_bsse[molecule][cluster][motion][config][method][basis]['no_vmfc'] = DE_no_vmfc
                        two_body_bsse[molecule][cluster][motion][config][method][basis]['vmfc'] = DE_vmfc
                        two_body_bsse[molecule][cluster][motion][config][method][basis]['bsse'] = DE_no_vmfc - DE_vmfc
    
    return two_body_bsse