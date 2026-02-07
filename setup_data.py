import os
import pandas as pd
import random

# Ensure data folder exists
if not os.path.exists('data'):
    os.makedirs('data')

print("🧪 INITIALIZING BIOPREDICTOR DATABASE...")

# --- 1. DEFINE REAL DATA BLOCKS ---
categories = {
    "Cancer_ABL1": {
        "seq": "MLEICLKLVGCKSKKGLSSSSSCYLEEALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAPKRNKPTVYGVSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKQGVRGAVSTLLQAPELPTKTRT",
        "drugs": [
            ("Imatinib", "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"),
            ("Nilotinib", "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5")
        ]
    },
    "Pain_COX2": {
        "seq": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "drugs": [
            ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
            ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
        ]
    },
    "Diabetes_AMPK": {
        "seq": "MAEKQKHDGRVKIGHYILGDTLGVGTFGKVKVGKHELTGHKVAVKILNRQKIRSLDVVGKIRREIQNLKLFRHPHIIKLYQVISTPTDFFMVMEYVSGGELFDYICKHGRVEEMEARRLFQQILSAVDYCHRHMVVHRDLKPENVLLDAHMNAKIADFGLSNMMSDGEFLRTSCGSPNYAAPEVISGRLYAGPEVDIWSSGVILYALLCGTLPFDDEHVPTLFKKIRGGVFYIPEYLNRSVATLLMHMLQVDPLKRATIKDIREHEWFKQDLPKYLFPEDPSYSSTMIDDEALKEVCEKFECSEEEVLSCLYNRNHQDPLAVAYHLIIDNRRIMNEAKDFYLATSPPDSFLDDHHLTR",
        "drugs": [
            ("Metformin", "CN(C)C(=N)NC(=N)N"),
            ("Glipizide", "CC1=CC=C(C=C1)S(=O)(=O)NC(=O)NC2CCCCC2")
        ]
    },
    "Mental_SERT": {
        "seq": "METTPLNSQKVLSECKDREDCQENGVLQKGVPTTADRAEPSQISNGYSNGVFHTRHASNVGFAWTEVLALALCLVFSLIGTLIGFGHSNSVRTVYLSSMALLTLGLVIWCSSSTTYLRYIIPHASLGKSGKLWQTLSICIFLGCTAALQVACLFQNCLLAWSAFQYIMIIAL",
        "drugs": [
            ("Prozac", "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F"),
            ("Zoloft", "CNC1CCC(C2=CC=CC=C21)C3=CC(=C(C=C3)Cl)Cl")
        ]
    },
    "Antibiotic_PBP": {
        "seq": "MKLKNTLGVVIGSLVAASAMNAFAQAQSKDGKALVAKDNVVSQALSKTQAADVDTVLDSLKTNGKPVTADKLVDLLNSKDIDTNLKAMTGGVGTAVADGKTVDTKVKTGDNVTGSSTAVSVSGNTLTGKTEGADLTVKDATGKVVGKAKVIDTGVAKDLTVTGTKTLDGKTVAVN",
        "drugs": [
            ("Penicillin", "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"),
            ("Amoxicillin", "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C")
        ]
    }
}

drugs_data = []
fasta_lines = []
interaction_data = []

# --- 2. GENERATE 300+ ENTRIES ---
print("   Generating files...")
count = 1

for cat_name, data in categories.items():
    # Create 15 variants of the protein to simulate biological diversity
    for p_var in range(15): 
        p_id = f"{cat_name}_v{p_var}"
        fasta_lines.append(f">{p_id}")
        fasta_lines.append(data['seq'])

        # Create drugs and interactions
        for d_name, smiles in data['drugs']:
             for d_var in range(5): # 5 batches per drug
                d_id = f"{d_name}_{count}"
                drugs_data.append({"drug_id": d_id, "smiles": smiles})
                
                # POSITIVE (Binds)
                interaction_data.append({"drug_id": d_id, "protein_id": p_id, "label": 1})
                
                # NEGATIVE (Mismatch - Picking a random other category)
                other_cat = random.choice(list(categories.keys()))
                if other_cat != cat_name:
                    bad_pid = f"{other_cat}_v0"
                    interaction_data.append({"drug_id": d_id, "protein_id": bad_pid, "label": 0})
                
                count += 1

# --- 3. SAVE TO FILES ---
# Save Drugs
pd.DataFrame(drugs_data).to_csv("data/drugs.csv", index=False)
# Save Interactions
pd.DataFrame(interaction_data).to_csv("data/interactions.csv", index=False)
# Save Fasta
with open("data/proteins.fasta", "w") as f:
    f.write("\n".join(fasta_lines))

print("✅ SUCCESS! Created 'data/drugs.csv', 'data/proteins.fasta', 'data/interactions.csv'.")