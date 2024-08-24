from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, Descriptors
from rdkit.Chem import rdMolDescriptors as Descriptors3D
import pandas as pd
import os

def get_morgan_fingerprints(mol, radius=2, nBits=1024):
    """
    extract "Morgan Fingerprints" using the new MorganGenerator
    """
    morgan_generator = AllChem.GetMorganGenerator(radius=radius, fpSize=nBits)
    return morgan_generator.GetFingerprint(mol)

def get_maccs_keys(mol):
    """
    extract "MACCS Keys" -> 166 binary template
    """
    return MACCSkeys.GenMACCSKeys(mol)

def get_molecular_descriptors(mol, is_3d=False):
    """
    Extract molecule descriptors based on the 2D or 3D configuration.
    """
    # 2D descriptors
    descriptors = {
        'MaxAbsEStateIndex': Descriptors.MaxAbsEStateIndex(mol),
        'MaxEStateIndex': Descriptors.MaxEStateIndex(mol),
        'MinAbsEStateIndex': Descriptors.MinAbsEStateIndex(mol),
        'MinEStateIndex': Descriptors.MinEStateIndex(mol),
        'qed': Descriptors.qed(mol),
        'SPS': Descriptors.SPS(mol),
        'MolWt': Descriptors.MolWt(mol),
        'HeavyAtomMolWt': Descriptors.HeavyAtomMolWt(mol),
        'ExactMolWt': Descriptors.ExactMolWt(mol),
        'NumValenceElectrons': Descriptors.NumValenceElectrons(mol),
        'NumRadicalElectrons': Descriptors.NumRadicalElectrons(mol),
        'MaxPartialCharge': Descriptors.MaxPartialCharge(mol),
        'MinPartialCharge': Descriptors.MinPartialCharge(mol),
        'MaxAbsPartialCharge': Descriptors.MaxAbsPartialCharge(mol),
        'MinAbsPartialCharge': Descriptors.MinAbsPartialCharge(mol),
        # 'FpDensityMorgan1': Descriptors.FpDensityMorgan1(mol),  # deprecation waring occurs 
        # 'FpDensityMorgan2': Descriptors.FpDensityMorgan2(mol),  # deprecation waring occurs 
        # 'FpDensityMorgan3': Descriptors.FpDensityMorgan3(mol),  # deprecation waring occurs 
        'BCUT2D_MWHI': Descriptors.BCUT2D_MWHI(mol), 
        'BCUT2D_MWLOW': Descriptors.BCUT2D_MWLOW(mol),
        'BCUT2D_CHGHI': Descriptors.BCUT2D_CHGHI(mol),
        'BCUT2D_CHGLO': Descriptors.BCUT2D_CHGLO(mol),
        'BCUT2D_LOGPHI': Descriptors.BCUT2D_LOGPHI(mol),
        'BCUT2D_LOGPLOW': Descriptors.BCUT2D_LOGPLOW(mol),
        'BCUT2D_MRHI': Descriptors.BCUT2D_MRHI(mol),
        'BCUT2D_MRLOW': Descriptors.BCUT2D_MRLOW(mol),
        'AvgIpc': Descriptors.AvgIpc(mol),
        'BalabanJ': Descriptors.BalabanJ(mol),
        'BertzCT': Descriptors.BertzCT(mol),
        'Chi0': Descriptors.Chi0(mol),
        'Chi0n': Descriptors.Chi0n(mol),
        'Chi0v': Descriptors.Chi0v(mol),
        'Chi1': Descriptors.Chi1(mol),
        'Chi1n': Descriptors.Chi1n(mol),
        'Chi1v': Descriptors.Chi1v(mol),
        'Chi2n': Descriptors.Chi2n(mol),
        'Chi2v': Descriptors.Chi2v(mol),
        'Chi3n': Descriptors.Chi3n(mol),
        'Chi3v': Descriptors.Chi3v(mol),
        'Chi4n': Descriptors.Chi4n(mol),
        'Chi4v': Descriptors.Chi4v(mol),
        'HallKierAlpha': Descriptors.HallKierAlpha(mol),
        'Ipc': Descriptors.Ipc(mol),
        'Kappa1': Descriptors.Kappa1(mol),
        'Kappa2': Descriptors.Kappa2(mol),
        'Kappa3': Descriptors.Kappa3(mol),
        'LabuteASA': Descriptors.LabuteASA(mol),
        'PEOE_VSA1': Descriptors.PEOE_VSA1(mol),
        'PEOE_VSA10': Descriptors.PEOE_VSA10(mol),
        'PEOE_VSA11': Descriptors.PEOE_VSA11(mol),
        'PEOE_VSA12': Descriptors.PEOE_VSA12(mol),
        'PEOE_VSA13': Descriptors.PEOE_VSA13(mol),
        'PEOE_VSA14': Descriptors.PEOE_VSA14(mol),
        'PEOE_VSA2': Descriptors.PEOE_VSA2(mol),
        'PEOE_VSA3': Descriptors.PEOE_VSA3(mol),
        'PEOE_VSA4': Descriptors.PEOE_VSA4(mol),
        'PEOE_VSA5': Descriptors.PEOE_VSA5(mol),
        'PEOE_VSA6': Descriptors.PEOE_VSA6(mol),
        'PEOE_VSA7': Descriptors.PEOE_VSA7(mol),
        'PEOE_VSA8': Descriptors.PEOE_VSA8(mol),
        'PEOE_VSA9': Descriptors.PEOE_VSA9(mol),
        'SMR_VSA1': Descriptors.SMR_VSA1(mol),
        'SMR_VSA10': Descriptors.SMR_VSA10(mol),
        'SMR_VSA2': Descriptors.SMR_VSA2(mol),
        'SMR_VSA3': Descriptors.SMR_VSA3(mol),
        'SMR_VSA4': Descriptors.SMR_VSA4(mol),
        'SMR_VSA5': Descriptors.SMR_VSA5(mol),
        'SMR_VSA6': Descriptors.SMR_VSA6(mol),
        'SMR_VSA7': Descriptors.SMR_VSA7(mol),
        'SMR_VSA8': Descriptors.SMR_VSA8(mol),
        'SMR_VSA9': Descriptors.SMR_VSA9(mol),
        'SlogP_VSA1': Descriptors.SlogP_VSA1(mol),
        'SlogP_VSA10': Descriptors.SlogP_VSA10(mol),
        'SlogP_VSA11': Descriptors.SlogP_VSA11(mol),
        'SlogP_VSA12': Descriptors.SlogP_VSA12(mol),
        'SlogP_VSA2': Descriptors.SlogP_VSA2(mol),
        'SlogP_VSA3': Descriptors.SlogP_VSA3(mol),
        'SlogP_VSA4': Descriptors.SlogP_VSA4(mol),
        'SlogP_VSA5': Descriptors.SlogP_VSA5(mol),
        'SlogP_VSA6': Descriptors.SlogP_VSA6(mol),
        'SlogP_VSA7': Descriptors.SlogP_VSA7(mol),
        'SlogP_VSA8': Descriptors.SlogP_VSA8(mol),
        'SlogP_VSA9': Descriptors.SlogP_VSA9(mol),
        'TPSA': Descriptors.TPSA(mol),
        'EState_VSA1': Descriptors.EState_VSA1(mol),
        'EState_VSA10': Descriptors.EState_VSA10(mol),
        'EState_VSA11': Descriptors.EState_VSA11(mol),
        'EState_VSA2': Descriptors.EState_VSA2(mol),
        'EState_VSA3': Descriptors.EState_VSA3(mol),
        'EState_VSA4': Descriptors.EState_VSA4(mol),
        'EState_VSA5': Descriptors.EState_VSA5(mol),
        'EState_VSA6': Descriptors.EState_VSA6(mol),
        'EState_VSA7': Descriptors.EState_VSA7(mol),
        'EState_VSA8': Descriptors.EState_VSA8(mol),
        'EState_VSA9': Descriptors.EState_VSA9(mol),
        'VSA_EState1': Descriptors.VSA_EState1(mol),
        'VSA_EState10': Descriptors.VSA_EState10(mol),
        'VSA_EState2': Descriptors.VSA_EState2(mol),
        'VSA_EState3': Descriptors.VSA_EState3(mol),
        'VSA_EState4': Descriptors.VSA_EState4(mol),
        'VSA_EState5': Descriptors.VSA_EState5(mol),
        'VSA_EState6': Descriptors.VSA_EState6(mol),
        'VSA_EState7': Descriptors.VSA_EState7(mol),
        'VSA_EState8': Descriptors.VSA_EState8(mol),
        'VSA_EState9': Descriptors.VSA_EState9(mol),
        'FractionCSP3': Descriptors.FractionCSP3(mol),
        'HeavyAtomCount': Descriptors.HeavyAtomCount(mol),
        'NHOHCount': Descriptors.NHOHCount(mol),
        'NOCount': Descriptors.NOCount(mol),
        'NumAliphaticCarbocycles': Descriptors.NumAliphaticCarbocycles(mol),
        'NumAliphaticHeterocycles': Descriptors.NumAliphaticHeterocycles(mol),
        'NumAliphaticRings': Descriptors.NumAliphaticRings(mol),
        'NumAromaticCarbocycles': Descriptors.NumAromaticCarbocycles(mol),
        'NumAromaticHeterocycles': Descriptors.NumAromaticHeterocycles(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHeteroatoms': Descriptors.NumHeteroatoms(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'NumSaturatedCarbocycles': Descriptors.NumSaturatedCarbocycles(mol),
        'NumSaturatedHeterocycles': Descriptors.NumSaturatedHeterocycles(mol),
        'NumSaturatedRings': Descriptors.NumSaturatedRings(mol),
        'RingCount': Descriptors.RingCount(mol),
        'MolLogP': Descriptors.MolLogP(mol),
        'MolMR': Descriptors.MolMR(mol),
        'fr_Al_COO': Descriptors.fr_Al_COO(mol),
        'fr_Al_OH': Descriptors.fr_Al_OH(mol),
        'fr_Al_OH_noTert': Descriptors.fr_Al_OH_noTert(mol),
        'fr_ArN': Descriptors.fr_ArN(mol),
        'fr_Ar_COO': Descriptors.fr_Ar_COO(mol),
        'fr_Ar_N': Descriptors.fr_Ar_N(mol),
        'fr_Ar_NH': Descriptors.fr_Ar_NH(mol),
        'fr_Ar_OH': Descriptors.fr_Ar_OH(mol),
        'fr_COO': Descriptors.fr_COO(mol),
        'fr_COO2': Descriptors.fr_COO2(mol),
        'fr_C_O': Descriptors.fr_C_O(mol),
        'fr_C_O_noCOO': Descriptors.fr_C_O_noCOO(mol),
        'fr_C_S': Descriptors.fr_C_S(mol),
        'fr_HOCCN': Descriptors.fr_HOCCN(mol),
        'fr_Imine': Descriptors.fr_Imine(mol),
        'fr_NH0': Descriptors.fr_NH0(mol),
        'fr_NH1': Descriptors.fr_NH1(mol),
        'fr_NH2': Descriptors.fr_NH2(mol),
        'fr_N_O': Descriptors.fr_N_O(mol),
        'fr_Ndealkylation1': Descriptors.fr_Ndealkylation1(mol),
        'fr_Ndealkylation2': Descriptors.fr_Ndealkylation2(mol),
        'fr_Nhpyrrole': Descriptors.fr_Nhpyrrole(mol),
        'fr_SH': Descriptors.fr_SH(mol),
        'fr_aldehyde': Descriptors.fr_aldehyde(mol),
        'fr_alkyl_carbamate': Descriptors.fr_alkyl_carbamate(mol),
        'fr_alkyl_halide': Descriptors.fr_alkyl_halide(mol),
        'fr_allylic_oxid': Descriptors.fr_allylic_oxid(mol),
        'fr_amide': Descriptors.fr_amide(mol),
        'fr_amidine': Descriptors.fr_amidine(mol),
        'fr_aniline': Descriptors.fr_aniline(mol),
        'fr_aryl_methyl': Descriptors.fr_aryl_methyl(mol),
        'fr_azide': Descriptors.fr_azide(mol),
        'fr_azo': Descriptors.fr_azo(mol),
        'fr_barbitur': Descriptors.fr_barbitur(mol),
        'fr_benzene': Descriptors.fr_benzene(mol),
        'fr_benzodiazepine': Descriptors.fr_benzodiazepine(mol),
        'fr_bicyclic': Descriptors.fr_bicyclic(mol),
        'fr_diazo': Descriptors.fr_diazo(mol),
        'fr_dihydropyridine': Descriptors.fr_dihydropyridine(mol),
        'fr_epoxide': Descriptors.fr_epoxide(mol),
        'fr_ester': Descriptors.fr_ester(mol),
        'fr_ether': Descriptors.fr_ether(mol),
        'fr_furan': Descriptors.fr_furan(mol),
        'fr_guanido': Descriptors.fr_guanido(mol),
        'fr_halogen': Descriptors.fr_halogen(mol),
        'fr_hdrzine': Descriptors.fr_hdrzine(mol),
        'fr_hdrzone': Descriptors.fr_hdrzone(mol),
        'fr_imidazole': Descriptors.fr_imidazole(mol),
        'fr_imide': Descriptors.fr_imide(mol),
        'fr_isocyan': Descriptors.fr_isocyan(mol),
        'fr_isothiocyan': Descriptors.fr_isothiocyan(mol),
        'fr_ketone': Descriptors.fr_ketone(mol),
        'fr_ketone_Topliss': Descriptors.fr_ketone_Topliss(mol),
        'fr_lactam': Descriptors.fr_lactam(mol),
        'fr_lactone': Descriptors.fr_lactone(mol),
        'fr_methoxy': Descriptors.fr_methoxy(mol),
        'fr_morpholine': Descriptors.fr_morpholine(mol),
        'fr_nitrile': Descriptors.fr_nitrile(mol),
        'fr_nitro': Descriptors.fr_nitro(mol),
        'fr_nitro_arom': Descriptors.fr_nitro_arom(mol),
        'fr_nitro_arom_nonortho': Descriptors.fr_nitro_arom_nonortho(mol),
        'fr_nitroso': Descriptors.fr_nitroso(mol),
        'fr_oxazole': Descriptors.fr_oxazole(mol),
        'fr_oxime': Descriptors.fr_oxime(mol),
        'fr_para_hydroxylation': Descriptors.fr_para_hydroxylation(mol),
        'fr_phenol': Descriptors.fr_phenol(mol),
        'fr_phenol_noOrthoHbond': Descriptors.fr_phenol_noOrthoHbond(mol),
        'fr_phos_acid': Descriptors.fr_phos_acid(mol),
        'fr_phos_ester': Descriptors.fr_phos_ester(mol),
        'fr_piperdine': Descriptors.fr_piperdine(mol),
        'fr_piperzine': Descriptors.fr_piperzine(mol),
        'fr_priamide': Descriptors.fr_priamide(mol),
        'fr_prisulfonamd': Descriptors.fr_prisulfonamd(mol),
        'fr_pyridine': Descriptors.fr_pyridine(mol),
        'fr_quatN': Descriptors.fr_quatN(mol),
        'fr_sulfide': Descriptors.fr_sulfide(mol),
        'fr_sulfonamd': Descriptors.fr_sulfonamd(mol),
        'fr_sulfone': Descriptors.fr_sulfone(mol),
        'fr_term_acetylene': Descriptors.fr_term_acetylene(mol),
        'fr_tetrazole': Descriptors.fr_tetrazole(mol),
        'fr_thiazole': Descriptors.fr_thiazole(mol),
        'fr_thiocyan': Descriptors.fr_thiocyan(mol),
        'fr_thiophene': Descriptors.fr_thiophene(mol),
        'fr_unbrch_alkane': Descriptors.fr_unbrch_alkane(mol),
        'fr_urea': Descriptors.fr_urea(mol),
    }

    if is_3d:
        # Ensure the molecule has explicit hydrogens
        mol = Chem.AddHs(mol)
        
        # Ensure the molecule has a 3D conformer
        if mol.GetNumConformers() == 0:
            try:
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.UFFOptimizeMolecule(mol)
            except:
                print(f"Failed to generate 3D conformer for molecule: {id}")
                return descriptors  # Return empty descriptors if 3D generation fails
        
        # 3D descriptors
        try:
            descriptors.update({
                '3D_CalcAUTOCORR3D': Descriptors3D.CalcAUTOCORR3D(mol),
                '3D_CalcAsphericity': Descriptors3D.CalcAsphericity(mol),
                '3D_CalcChi0n': Descriptors3D.CalcChi0n(mol),
                '3D_CalcChi0v': Descriptors3D.CalcChi0v(mol),
                '3D_CalcChi1n': Descriptors3D.CalcChi1n(mol),
                '3D_CalcChi1v': Descriptors3D.CalcChi1v(mol),
                '3D_CalcChi2n': Descriptors3D.CalcChi2n(mol),
                '3D_CalcChi2v': Descriptors3D.CalcChi2v(mol),
                '3D_CalcChi3n': Descriptors3D.CalcChi3n(mol),
                '3D_CalcChi3v': Descriptors3D.CalcChi3v(mol),
                '3D_CalcChi4n': Descriptors3D.CalcChi4n(mol),
                '3D_CalcChi4v': Descriptors3D.CalcChi4v(mol),
                '3D_CalcCoulombMat': Descriptors3D.CalcCoulombMat(mol),
                '3D_CalcCrippenDescriptors': Descriptors3D.CalcCrippenDescriptors(mol)[0],
                '3D_CalcEEMcharges': Descriptors3D.CalcEEMcharges(mol),
                '3D_CalcEccentricity': Descriptors3D.CalcEccentricity(mol),
                '3D_CalcExactMolWt': Descriptors3D.CalcExactMolWt(mol),
                '3D_CalcFractionCSP3': Descriptors3D.CalcFractionCSP3(mol),
                # '3D_CalcGETAWAY': Descriptors3D.CalcGETAWAY(mol),  # stuck in some SMILES code (can't catch pattern) 
                '3D_CalcHallKierAlpha': Descriptors3D.CalcHallKierAlpha(mol),
                '3D_CalcInertialShapeFactor': Descriptors3D.CalcInertialShapeFactor(mol),
                '3D_CalcKappa1': Descriptors3D.CalcKappa1(mol),
                '3D_CalcKappa2': Descriptors3D.CalcKappa2(mol),
                '3D_CalcKappa3': Descriptors3D.CalcKappa3(mol),
                '3D_CalcLabuteASA': Descriptors3D.CalcLabuteASA(mol),
                '3D_CalcMORSE': Descriptors3D.CalcMORSE(mol),
                '3D_CalcMolFormula': Descriptors3D.CalcMolFormula(mol),
                '3D_CalcNPR1': Descriptors3D.CalcNPR1(mol),
                '3D_CalcNPR2': Descriptors3D.CalcNPR2(mol),
                '3D_CalcNumAliphaticCarbocycles': Descriptors3D.CalcNumAliphaticCarbocycles(mol),
                '3D_CalcNumAliphaticHeterocycles': Descriptors3D.CalcNumAliphaticHeterocycles(mol),
                '3D_CalcNumAliphaticRings': Descriptors3D.CalcNumAliphaticRings(mol),
                '3D_CalcNumAmideBonds': Descriptors3D.CalcNumAmideBonds(mol),
                '3D_CalcNumAromaticCarbocycles': Descriptors3D.CalcNumAromaticCarbocycles(mol),
                '3D_CalcNumAromaticHeterocycles': Descriptors3D.CalcNumAromaticHeterocycles(mol),
                '3D_CalcNumAromaticRings': Descriptors3D.CalcNumAromaticRings(mol),
                # '3D_CalcNumAtomStereoCenters': Descriptors3D.CalcNumAtomStereoCenters(mol),  # err: numUnspecifiedStereoCenters called without stereo being assigned
                '3D_CalcNumAtoms': Descriptors3D.CalcNumAtoms(mol),
                '3D_CalcNumBridgeheadAtoms': Descriptors3D.CalcNumBridgeheadAtoms(mol),
                '3D_CalcNumHBA': Descriptors3D.CalcNumHBA(mol),
                '3D_CalcNumHBD': Descriptors3D.CalcNumHBD(mol),
                '3D_CalcNumHeavyAtoms': Descriptors3D.CalcNumHeavyAtoms(mol),
                '3D_CalcNumHeteroatoms': Descriptors3D.CalcNumHeteroatoms(mol),
                '3D_CalcNumHeterocycles': Descriptors3D.CalcNumHeterocycles(mol),
                '3D_CalcNumLipinskiHBA': Descriptors3D.CalcNumLipinskiHBA(mol),
                '3D_CalcNumLipinskiHBD': Descriptors3D.CalcNumLipinskiHBD(mol),
                '3D_CalcNumRings': Descriptors3D.CalcNumRings(mol),
                '3D_CalcNumRotatableBonds': Descriptors3D.CalcNumRotatableBonds(mol),
                '3D_CalcNumSaturatedCarbocycles': Descriptors3D.CalcNumSaturatedCarbocycles(mol),
                '3D_CalcNumSaturatedHeterocycles': Descriptors3D.CalcNumSaturatedHeterocycles(mol),
                '3D_CalcNumSaturatedRings': Descriptors3D.CalcNumSaturatedRings(mol),
                '3D_CalcNumSpiroAtoms': Descriptors3D.CalcNumSpiroAtoms(mol),
                # '3D_CalcNumUnspecifiedAtomStereoCenters': Descriptors3D.CalcNumUnspecifiedAtomStereoCenters(mol),  # err: numUnspecifiedStereoCenters called without stereo being assigned
                '3D_CalcOxidationNumbers': Descriptors3D.CalcOxidationNumbers(mol),
                '3D_CalcPBF': Descriptors3D.CalcPBF(mol),
                '3D_CalcPMI1': Descriptors3D.CalcPMI1(mol),
                '3D_CalcPMI2': Descriptors3D.CalcPMI2(mol),
                '3D_CalcPMI3': Descriptors3D.CalcPMI3(mol),
                '3D_CalcPhi': Descriptors3D.CalcPhi(mol),
                '3D_CalcRDF': Descriptors3D.CalcRDF(mol),
                '3D_CalcRadiusOfGyration': Descriptors3D.CalcRadiusOfGyration(mol),
                '3D_CalcSpherocityIndex': Descriptors3D.CalcSpherocityIndex(mol),
                '3D_CalcTPSA': Descriptors3D.CalcTPSA(mol),
                '3D_CalcWHIM': Descriptors3D.CalcWHIM(mol)
            })
        except Exception as e:
            print(f"Failed to calculate 3D descriptors for molecule: {id}, error: {e}")

    return descriptors

def smiles_to_features(rawdata_filepath, save_dirpath):    
    # init
    morgan_data = []  # for Morgan fingerprints DataFrame
    maccs_data = []   # for MACCS keys DataFrame
    descriptor_data = []  # for molecular descriptors DataFrame

    df = pd.read_csv(rawdata_filepath)
    row_count = df.shape[0]
    counter = 1
    for smi, id in zip(df["Smiles"], df["Molecule ChEMBL ID"]):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            print(f"{counter}/{row_count}: {id}")

            # Morgan fingerprints
            morgan_fp = get_morgan_fingerprints(mol)
            morgan_fp_list = list(morgan_fp)
            morgan_data.append(morgan_fp_list)

            # MACCS keys
            maccs_keys = get_maccs_keys(mol)
            maccs_keys_list = list(maccs_keys)[1:]  # del 0 index value
            maccs_data.append(maccs_keys_list)
            
            # Molecular descriptors
            descriptors = get_molecular_descriptors(mol, is_3d=True)  # contain 3d 
            descriptor_values = list(descriptors.values())
            descriptor_data.append(descriptor_values)
    
            counter+=1
    
    # Generate column names for each DataFrame
    morgan_columns = [f'Morgan_{i}' for i in range(1024)]
    maccs_columns = [f'MACCS_{i}' for i in range(1, 167)]
    descriptor_columns = list(descriptors.keys())

    # Create and save DataFrames for each set of features
    morgan_df = pd.DataFrame(morgan_data, columns=morgan_columns)
    morgan_df = pd.concat([df.copy(), morgan_df], axis=1)
    morgan_df.to_csv(os.path.join(save_dirpath, "train_rdkit_morgan.csv"), index=False)
    print(f"Morgan fingerprints saved to {os.path.join(save_dirpath, 'train_rdkit_morgan.csv')}")

    maccs_df = pd.DataFrame(maccs_data, columns=maccs_columns)
    maccs_df = pd.concat([df.copy(), maccs_df], axis=1)
    maccs_df.to_csv(os.path.join(save_dirpath, "train_rdkit_maccs.csv"), index=False)
    print(f"MACCS keys saved to {os.path.join(save_dirpath, 'train_rdkit_maccs.csv')}")

    descriptor_df = pd.DataFrame(descriptor_data, columns=descriptor_columns)
    descriptor_df = pd.concat([df.copy(), descriptor_df], axis=1)
    descriptor_df.to_csv(os.path.join(save_dirpath, "train_rdkit_descriptors.csv"), index=False)
    print(f"Molecular descriptors saved to {os.path.join(save_dirpath, 'train_rdkit_descriptors.csv')}")

# execute
if __name__ == "__main__":
    RAWDATA_FILEPATH = r"data\preprocessed\train_extracted_SMILES.csv"
    SAVE_DIRPATH = r"data\preprocessed"
    smiles_to_features(RAWDATA_FILEPATH, SAVE_DIRPATH)