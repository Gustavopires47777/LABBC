from Bio import Entrez
from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
import matplotlib.pyplot as plt
from Bio.Phylo import PhyloXML
from Bio.Phylo.Consensus import *
import pylab # type: ignore
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyScorer, NNITreeSearcher,ParsimonyTreeConstructor

align2 = AlignIO.read('ExemploArquivo.phy','phylip') #Arquivo precisa estar alinhado para ser lido
print(align2, '\n') #printar alinhamento armazenado

#DISTANCIA GENÉTICA 
calculator = DistanceCalculator('identity') #Classe do módulo Phylo que armazena a distancia genética
dist_mat = calculator.get_distance(align2) #método get_distance para exibir a distancia genética da variável align2, a qual armazena o seu arquivo de alinhamento
print (dist_mat) #printar distancia de alinhamento

#CONSTRUÇÃO DE ARVORE FILOGENÉTICA COM O ARQUIVO DE ALINHAMENTO
constructor = DistanceTreeConstructor(calculator, 'nj') #metodo para construção de arvore calculando a partir da distancia genética entre os alinhamentos # 'nj' / 'upgma' tree construction approaches 
tree_0 = constructor.build_tree(align2)
print(tree_0) #printar arvore filogenetica
#Phylo.draw(tree_0)

#CÁLCULO DE PARSIMONIA
#ParsimonyScorer calcula o score de parsimonia da árvore target através do input de alinhamento, e TreeSearcher é o metodo que pesquisa a melhor configuração de árvore
scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
# constructor = ParsimonyTreeConstructor(searcher, starting_tree)
pars_tree = constructor.build_tree(align2)
print(pars_tree)

#CONSENSUS TREE, modulo importado do Bio.Phylo.Consensus
#Pode evocar os métodos stric_consensus, majority_consensus e adam_consensus
trees = list(Phylo.parse("seu_arquivo.nwk","newick")) #exemplo de formato
strict_tree = strict_tree(trees)
majority_tree = majority_consensus(trees,0.5) #teste de curva ROC verdadeiro-falso - razao de verossimilhança
adam_tree = adam_consensus(trees)

'''Instead of using 50% as the cutoff, the majority_consensus method allows you to set your own one 
by providing an extra cutoff parameter(0~1, 0 by default). That means it can also work as strict consensus 
algorithm when the cutoff equals 1. One more thing different to strict and adam consensus tree, the result 
majority rule consensus tree has branch support value that are automatically assigned during calculation.'''

#Bootstrap Methods
#Para construir uma árvore consensus é necessário construir uma lista de replicatas de arvores, para isso se utiliza Bio.Phylo.Consensus
msa = AlignIO.read("seu_arquivo.phy","phylip")
msas = bootstrap(msa, 100)

#Aceita metodos de alinhamento múltiplo, quando é setado 100, a arvore é replicada 100 vezes

calculator = DistanceCalculator("blosum62")
constructor = DistanceTreeConstructor(calculator)
trees = bootstrap_trees(msa,100, constructor)

#Another useful function to replicate trees as a extra parameter called consensus_tree
consensus_tree = bootstrap_consensus(msa, 100, constructor, majority_consensus)

#CONVERSÃO DE FORMATO DE ARQUIVO
tree1 = Phylo.read("meu_exemplo.nwk", "newick")
Phylo.convert("meu_exemplo.nwk", "newick", "meu_exemplo.xml", "phyloxml") #Converter formato e extensão
