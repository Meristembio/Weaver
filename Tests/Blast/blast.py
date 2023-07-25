import os
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
from pyblast import BioBlast
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pyblast.utils import make_linear, make_circular

warnings.simplefilter('ignore', BiopythonWarning)

plasmids_path = "plasmids"
query_string = "aataataattttttcatgttgaaaatctccaaaaaaaaaggctccaaaaggagcctttaattgtatcggtttatcagcttgctttttatgacaacttgacggctacatcattcactttttcttcacaaccggcacggaactcgctcgggctggccccggtgcattttttaaatacccgcgagaaatagagttgatcgtcaaaaccaacattgcgaccgacggtggcgataggcatccgggtggtgctcaaaagcagcttcgcctggctgatacgttggtcctcgcgccagcttaagacgctaatccctaactgctggcggaaaagatgtgacagacgcgacggcgacaagcaaacatgctgtgcgacgctggcgatatcaaaattgctgtctgccaggtgatcgctgatgtactgacaagcctcgcgtacccgattatccatcggtggatggagcgactcgttaatcgcttccatgcgccgcagtaacaattgctcaagcagatttatcgccagcagctccgaatagcgcccttccccttgcccggcgttaatgatttgcccaaacaggtcgctgaaatgcggctggtgcgcttcatccgggcgaaagaaccccgtattggcaaatattgacggccagttaagccattcatgccagtaggcgcgcggacgaaagtaaacccactggtgataccattcgcgagcctccggatgacgaccgtagtgatgaatctctcctggcgggaacagcaaaatatcacccggtcggcaaacaaattctcgtccctgatttttcaccaccccctgaccgcgaatggtgagattgagaatataacctttcattcccagcggtcggtcgataaaaaaatcgagataaccgttggcctcaatcggcgttaaacccgccaccagatgggcattaaacgagtatcccggcagcaggggatcattttgcgcttcagccatacttttcatactcccgccattcagagaagaaaccaattgtccatattgcatcagacattgccgtcactgcgtcttttactggttcttctcgctaaccaaaccggtaaccccgcttattaaaagcattctgtaacaaagcgggaccaaagccatgacaaaaacgcgtaacaaaagtgtctataatcacggcagaaaagtccacattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgacgctttttatcgcaactctctactgtttctccatacccgtttttttgggctagc"

queries = [
  SeqRecord(Seq(query_string))
]
queries = make_linear(queries)

print("Creating subjects")
subjects = []
for file in os.listdir(plasmids_path):
    if file.endswith('.gb'):
        try:
            record = list(SeqIO.parse(os.path.join(plasmids_path, file), 'gb'))[0]
            record.id = file.__str__()
            subjects.append(record)
        except ValueError:
            pass
        except AttributeError:
            pass
        except IndexError:
            pass
    if len(subjects) > 1000:
        break

subjects = make_circular(subjects)

print("Making blast")
blast = BioBlast(subjects, queries)
results = blast.blastn()
print(results)

'''
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

query = record.seq.upper()

alignments = []
for i, plasmid in enumerate(plasmids):
    print("Processing plasmid #" + i.__str__())
    try:
        plasmid_seq = list(plasmid)[0].seq.upper()
        if plasmid_seq and query:
            if aligner.score(plasmid_seq, query) > score_threshold:
                alignment_result = aligner.align(plasmid_seq, query)
                if len(alignment_result):
                    optimal_alignment = next(alignment_result)
                    alignments.append({
                        'plasmid': plasmid,
                        'score': optimal_alignment.score,
                        'to_str': optimal_alignment
                    })
    except:
        continue
alignments = sorted(alignments, key=lambda d: d['score'], reverse=True)
print(alignments)
'''