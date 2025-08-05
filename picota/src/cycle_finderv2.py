import collections
from collections import deque
from src.cycle_kmer_hash import filter_cycles_with_kmer
from src.cycle_kmer_hash import print_progress_bar

class Graph(object):
    def __init__(self, edges):
        self.edges = edges
        self.adj = Graph._build_adjacency_list(edges)

    @staticmethod
    def _build_adjacency_list(edges):
        adj = collections.defaultdict(list)
        for edge in edges:
            adj[edge[0]].append(edge[1])
        return adj

class GraphWork:
    def __init__(self):
        self.cycles = []
        self.reverse_or_cycles = []
        self.paths = [] #no option reverse for paths
        self.allPaths = []
        self.reverseallPaths = []
        self.visited = set()
        self.find_all_path = False
        self.path_limit = 15

    def dfs(self, G):
        discovered = set()
        finished = set()
        
        for u in G.adj.copy():
            if u not in discovered and u not in finished:
                discovered, finished = self.dfs_visit(G, u, discovered, finished)

    def dfs_iterative(self, G):
        # Başlangıç durumu
        stack = []
        discovered = set()  # Ziyaret edilen düğümler
        finished = set()    # Tamamen işlenmiş düğümler
        self.cycles = []    # Tespit edilen döngüler
        self.reverse_or_cycles = []  # Ters döngüler

        # Grafın tüm düğümleri üzerinden iterasyon
        for start_node in G.adj.copy():
            if start_node not in discovered:
                stack.append((start_node, iter(G.adj[start_node])))

                while stack:
                    current_node, neighbors = stack[-1]  # Yığının tepesindeki düğüm

                    if current_node not in discovered:
                        discovered.add(current_node)

                    try:
                        neighbor = next(neighbors)  # Sıradaki komşuyu al

                        # Döngü kontrolü
                        if neighbor in discovered:
                            self.cycles.append((current_node, neighbor))
                            if self.find_all_path:
                                self.capturePaths(G, neighbor, current_node, reverse=False)
                            else:
                                self.find_paths(neighbor, current_node, G)
                            continue

                        elif neighbor.replace('-', '+') in discovered or neighbor.replace('+', '-') in discovered:
                            self.reverse_or_cycles.append((current_node, neighbor))
                            the_changed_v = (
                                neighbor.replace('-', '+') if neighbor.replace('-', '+') in discovered
                                else neighbor.replace('+', '-')
                            )
                            if self.find_all_path:
                                self.capturePaths(G, the_changed_v, current_node, reverse=True)
                            else:
                                self.find_paths(the_changed_v, current_node, G)
                            continue

                        if neighbor not in finished:
                            # Yığının üstüne yeni komşuyu ekle
                            stack.append((neighbor, iter(G.adj[neighbor])))

                    except StopIteration:
                        # Komşular bittiğinde yığından çıkar
                        stack.pop()
                        discovered.remove(current_node)
                        finished.add(current_node)

        return discovered, finished




    def dfs_visit(self, G, u, discovered, finished):
        discovered.add(u)
        
        for v in G.adj[u]:
            # Detect cycles

            if v in discovered:
                #print(f"Cycle detected: found a back edge from {u} to {v}.")
                self.cycles.append((u,v))
                #self.find_paths(v, u, G)
                #exponential 
                
                if self.find_all_path:
                    self.capturePaths(G, v, u, reverse=False)
                else:
                    self.find_paths(v, u, G)
                break
            elif v.replace('-','+') in discovered or v.replace('+','-') in discovered:
                self.reverse_or_cycles.append((u,v))

                if v.replace('-','+') in discovered:
                    the_changed_v = v.replace('-','+')
                else:
                    the_changed_v = v.replace('+','-')

                if self.find_all_path:
                    self.capturePaths(G, the_changed_v, u, reverse=True)
                else:
                    self.find_paths(the_changed_v, u, G)
                break

            # Recurse into DFS tree
            if v not in finished:
                self.dfs_visit(G, v, discovered, finished)

        discovered.remove(u)
        finished.add(u)

        return discovered, finished


    #Careful it is exponential time because it is NP Hard
    def findAllPosPaths(self, graph, src, dest, path, reverse):
        
        self.visited.add(src)
        path.append(src)

        if src == dest and reverse == False:
            self.allPaths.append(path.copy())
        elif src == dest and reverse == True:
            self.reverseallPaths.append(path.copy())

        elif len(path) >= self.path_limit:
            LIMIT = True
        else:
            for i in graph.adj[src]:
                if i not in self.visited:
                    self.findAllPosPaths(graph, i, dest, path, reverse)
        
        path.pop()
        self.visited.remove(src)

    def capturePaths(self, graph, src, dest, reverse):
        path = []
        self.visited = set()
        try:
            self.findAllPosPaths(graph, src, dest, path, reverse)
        except RecursionError:
            print('Recursion Error happened try with another settings, finding all paths can not be possible!')
        except Exception as e:
            print(e,'Unkown error happened, , finding all paths can not be possible!, try with simple path setting!')


    def isReachable(self, graph, src, dest, discovered, path):
    
        # mark the current node as discovered
        discovered.add(src)
    
        # include the current node in the path
        path.append(src)
    
        # if destination vertex is found
        if src == dest:
            return True
    
        # do for every edge (src, i)
        for i in graph.adj[src]:
    
            # if `u` is not yet discovered
            if i not in discovered:
                # return true if the destination is found
                if self.isReachable(graph, i, dest, discovered, path):
                    return True
    
        # backtrack: remove the current node from the path
        path.pop()
    
        # return false if destination vertex is not reachable from src
        return False


    def find_paths(self, src,dest,the_graph):
            # to keep track of whether a vertex is discovered or not
        discovered = set()
    
        # List to store the complete path between source and destination
        path = deque()
    
        # perform DFS traversal from the source vertex to check the connectivity
        # and store path from the source vertex to the destination vertex
        if self.isReachable(the_graph, src, dest, discovered, path):
            #print(f'Path exists from vertex {src} to vertex {dest}')
            #print(f'The complete path is', list(path))
            self.paths.append(list(path))



    def reverse_complement(self, seq):
        complement_dict = {"A":"T", "T":"A", "G":"C", "C":"G", "a":"t", "t":"a", "g":"c", "c":"g"}
        return "".join([complement_dict[nucleotide] for nucleotide in reversed(seq)])


    def reverse_sign(self, sign):
        sign_dict = {"+" : "-", "-" : "+"}
        return sign_dict[sign]


    def parse_gfa(self, path_to_gfa):
        edge_dict = {}
        node_dict = {}
        print(path_to_gfa)
        with open(path_to_gfa, "r") as gfa:

            for line in gfa:

                #Add sequences as nodes, one node for each strand of the DNA
                if line.startswith("S"):
                    temp = line.strip("\n").split("\t")
                    node_dict[temp[1] + "+"] = {"Name": temp[1] + "+", "Sequence": temp[2]}
                    node_dict[temp[1] + "-"] = {"Name": temp[1] + "-", "Sequence": self.reverse_complement(temp[2])}

                # Add links as edges, one edge for each link. Since two nodes for each sequence exits reverse links are also
                # generated and added as edges.
                elif line.startswith("L"):
                    temp = line.strip("\n").split("\t")
                    if (temp[1] + temp[2], temp[3] + temp[4]) not in edge_dict.keys() and (
                            temp[3] + self.reverse_sign(temp[4]), temp[1] + self.reverse_sign(temp[2])) not in edge_dict.keys():
                    
                        edge_dict[(temp[1] + temp[2], temp[3] + temp[4])] = {"From": temp[1] + temp[2],
                                                                            "To": temp[3] + temp[4], "FromOrient": temp[2],
                                                                            "ToOrient": temp[4], "Overlap": temp[5]}

                        edge_dict[(temp[3] + self.reverse_sign(temp[4]), temp[1] + self.reverse_sign(temp[2]))] = {
                            "From": temp[3] + self.reverse_sign(temp[4]), "To": temp[1] + self.reverse_sign(temp[2]),
                            "FromOrient": self.reverse_sign(temp[4]), "ToOrient": self.reverse_sign(temp[2]), "Overlap": temp[5]}
        return node_dict, edge_dict


    def generate_genome_graph(self, node_dict, edge_dict):
        

        #add nodes and edges to the directed graph
        #G.add_nodes_from(list(node_dict.keys()))
        #nx.set_node_attributes(G, node_dict)

        #G.add_edges_from(list(edge_dict.keys()))
        #nx.set_edge_attributes(G, edge_dict)

        G = Graph(list(edge_dict.keys()))

        return G


class Cycle:
    def __init__(self, name, sequence, length, component_number, path):
        self.name = name
        self.sequence = sequence
        self.length = length
        self.component_number = component_number
        self.path = path
        self.reverseOriented = False


def reverse_complement(seq):
    complement_dict = {"A":"T", "T":"A", "G":"C", "C":"G", "a":"t", "t":"a", "g":"c", "c":"g"}
    return "".join([complement_dict[nucleotide] for nucleotide in reversed(seq)])


def find_overlap_length(overlap):
    if overlap == "0M":
        overlap_length = 0
    else:
        overlap_length = int(overlap[:-1])
    return overlap_length



def is_similar_polynomial(seq1, seq2, k_mer_sim, threshold_sim):
    # Precompute k-mers from seq1 and its reverse complement
    kmer_set = set(
        seq1[i:i+k_mer_sim] for i in range(len(seq1) - k_mer_sim + 1)
    ) | set(
        reverse_complement(seq1)[i:i+k_mer_sim] for i in range(len(seq1) - k_mer_sim + 1)
    )

    # Sliding window to count matching k-mers in seq2
    exact_count = 0
    for i in range(len(seq2) - k_mer_sim + 1):
        kmer = seq2[i:i+k_mer_sim]
        if kmer in kmer_set:
            exact_count += 1

    # Check threshold similarity
    total_kmers = len(seq2) - k_mer_sim + 1
    if (exact_count / total_kmers) * 100 > threshold_sim:
        return True
    return False


#%95 threshold,
def is_similar(seq1, seq2, k_mer_sim, threshold_sim):
    C_t=list(map(''.join, zip(*[iter(seq2)]*k_mer_sim)))
    exact_count = 0
    for el in C_t:
        if el in seq1 or el in reverse_complement(seq1):
            exact_count += 1
    if (exact_count/len(C_t))*100 > threshold_sim:
        return True
    return False


#path = [340+, 350-, 20-]
#nodes_len = {340+: 245}
#new_parts = [[340+, 350-, 20-],[340+, 350-, 10-]]
def cycle_match_based_on_contig_id(path, nodes_len, new_parts, threshold=70):
    # İlk path için zaten new_parts boş, bu yüzden direkt olarak ekleyebiliriz.
    if not new_parts:
        return True

    # Path'i new_parts içinde yer alan her path ile karşılaştır ve benzerliği kontrol et
    path_set = set(path)
    #path_length = sum(nodes_len.get(node, 0) for node in path)

    
    path_length = sum(nodes_len.get(node, 0) for node in path)

    for part in new_parts:
    
        existing_path_set = set(part)
        common_nodes = path_set.intersection(existing_path_set)

        common_length = sum(nodes_len[node] for node in common_nodes if node in nodes_len)
        part_length = sum(nodes_len.get(node, 0) for node in part)

        max_length = max(part_length, path_length)
        
        if max_length > 0:
            identity = common_length / max_length * 100
        else:
            identity = 0  
        
        if identity > threshold:
            #print(identity)
            #print(path)
            #print(existing_path_set)
            #print(common_length, max_length)
            #input()
            return False  # Eğer benzerlik eşik değerini aşarsa, bu path'i eklemeyiz.

    return True


def cycle_info_optimized(path, nodes, edges, cycle_info_list):   
    component_number = len(path)
    total_lem = 0
    the_final_seq = ''

    # Prepare sequences and overlap checks upfront
    sequences = [nodes[path[i]]["Sequence"] for i in range(component_number)]
    overlaps = [edges[path[i], path[i + 1]]["Overlap"] for i in range(component_number - 1)]

    # Build the final sequence
    for i in range(component_number - 1):
        overlap = overlaps[i]
        if overlap == "*":
            print(f"At least one of the edges has a non-specified overlap in this cycle: {path}\nSkipping cycle...")
            return None
    
        overlap_length = find_overlap_length(overlap)
        total_lem += len(sequences[i])
        
        # Append sequence with or without overlap
        the_final_seq += sequences[i][:-overlap_length] if overlap_length else sequences[i]

    # Append the last sequence after processing the penultimate
    the_final_seq += sequences[-1]

    # Reverse complement calculation done only once
    reverse_final_seq = reverse_complement(the_final_seq)

    # Check for exact match with existing cycles
    for elm in cycle_info_list:
        if (elm.sequence == the_final_seq or elm.sequence == reverse_final_seq):
            return 'Pass'
    
    # Return a new Cycle object if no match is found
    return Cycle('name', the_final_seq, len(the_final_seq), component_number, path)


def cycle_info(path, nodes, edges, cycle_info_list):

    component_number = len(path)
    
    the_final_seq = ''
    total_lem = 0

    #Self Circuit İdentification
    if component_number == 1:
        the_final_seq = nodes[path[0]]["Sequence"]
        the_edge = edges[path[0],path[0]]
        overlap = the_edge["Overlap"]
        if overlap == "*":
            print(f"At least on of the edges has non-specified overlap in this cycle: {cycle}\n Skipping cycle...")
            return None
        overlap_length = find_overlap_length(overlap)
    else:
        for iterator in range (0, component_number-1):
            the_edge = edges[path[iterator],path[iterator+1]]
            first_node_seq = nodes[path[iterator]]["Sequence"]
            second_node_seq = nodes[path[iterator+1]]["Sequence"]
            total_lem += len(first_node_seq)
            overlap = the_edge["Overlap"]
            if overlap == "*":
                print(f"At least on of the edges has non-specified overlap in this cycle: {cycle}\n Skipping cycle...")
                return None
            overlap_length = find_overlap_length(overlap)
            if overlap_length != 0:
                the_final_seq += first_node_seq[:-overlap_length]
            else:
                the_final_seq += first_node_seq
            if iterator == component_number-2:
                the_final_seq += second_node_seq

    for elm in cycle_info_list:


       if overlap_length != 0:
           if elm.sequence in (the_final_seq[:-overlap_length] + the_final_seq) or elm.sequence in (reverse_complement(the_final_seq)[:-overlap_length] + reverse_complement(the_final_seq)):
               #There is an exact cycle so dont take it
               return 'Pass'
       else:
           if elm.sequence in (the_final_seq + the_final_seq) or elm.sequence in (reverse_complement(the_final_seq) + reverse_complement(the_final_seq)):
               #There is an exact cycle so dont take it
               return 'Pass'


    return Cycle('name', the_final_seq, len(the_final_seq), component_number, path)






def cycle_analysis(path_to_data, out_cycle_file, find_all_path, path_limit, min_size_of_cycle, max_size_of_cycle, name_prefix_cycle, min_component_number, max_component_number, k_mer_sim, threshold_sim):

    
    #path_to_data = "test_data/assembly/gfa_files/test_data_P.nitroreducens_1Iinsides.gfa"
    #path_to_data = "test_data/k39.gfa"
    #path_to_data = "CompTestData/TnXc5-k119.gfa"
    #path_to_data = "sra_download/assembly/gfa_files/k21.gfa"


    #find_all_path = True
    #path_limit = 25

    #out_cycle_file = "test_data/cycles/cycles79Ec.fasta"

    #min_size_of_cycle = 3000
    #max_size_of_cycle = 100000
    #name_prefix_cycle = 'Cycle'
    #min_component_number = 1
    #max_component_number = 25


    #k_mer_sim = 200
    #threshold_sim = 99



    ###########################################
    #############################################


    GW = GraphWork()
    GW.find_all_path = find_all_path
    GW.path_limit = path_limit

    print('Parsing the GFA File')
    node_dict, edge_dict = GW.parse_gfa(path_to_data)
    print('Generating the Genome Graph')
    genome_graph = GW.generate_genome_graph(node_dict, edge_dict)
    print('Detecting the Cycles')
    try:
        #GW.dfs(genome_graph) it is not iterative can ve problematic
        GW.dfs_iterative(genome_graph)
    except RecursionError:
        print('ERROR(Recursion): DFS can not be possible with this graph, try to assemble with higher k-mer!')        
        f = open(out_cycle_file, "w")
        f.write('Error recursion')
        f.close()
        raise RecursionError('input is too big')
        return 'ERROR'
    except Exception as e:
        f = open(out_cycle_file, "w")
        f.write('Error recursion')
        f.close()
        print('ERROR: Unkown error happened, finding cycles maybe are not possible\n', e)
        raise RecursionError('input is too big')
        return 'ERROR'


    print('Finding Paths from cycles...')
    print('Cycle number: ',len(GW.cycles))
    print('All Paths: ', len(GW.allPaths))
    print('Reverse Paths: ', len(GW.reverseallPaths))

    name_it = 0
    cycle_info_list = []

    #selected paths
    if find_all_path:
        selected_paths = GW.allPaths
    else:
        selected_paths = GW.paths
    
    print('Len of Cycles before elimination: ', len(selected_paths))
    i_count = 0
    total_paths = len(selected_paths)
    
    new_paths = []

    # node_dict örneğinizin içerdiği yapıya uygun olarak
    node_lengths = {key: len(value["Sequence"]) for key, value in node_dict.items()}

    print('Contig ID based -')
    print('Before Elimination:', len(selected_paths))
    p_count = 0
    for path in selected_paths:
        p_count += 1
        print_progress_bar(p_count, len(selected_paths), prefix='Processing:', suffix='Complete')
        if cycle_match_based_on_contig_id(path, node_lengths, new_paths):
            new_paths.append(path)
    

    for path in new_paths:
        i_count += 1

        print_progress_bar(i_count, len(new_paths), prefix='Processing:', suffix='Complete')

        
        cycle_inf_obj = cycle_info_optimized(path, node_dict, edge_dict, cycle_info_list)
        if cycle_inf_obj == None:
            print('An error occured')
            continue
        
        if cycle_inf_obj == 'Pass':
            #print('Exact Cycle is existed')
            continue
        

        if cycle_inf_obj.length >= min_size_of_cycle and cycle_inf_obj.length <= max_size_of_cycle and cycle_inf_obj.component_number >= min_component_number and cycle_inf_obj.component_number <= max_component_number:
            name_it += 1
            cycle_inf_obj.name = name_prefix_cycle + '_' + str(name_it)
            cycle_info_list.append(cycle_inf_obj)

    print('\nAfter Elimination',len(cycle_info_list))
    print(cycle_info_list)
    print('\nReverse Paths reducing:')
    i_count_2 = 0
    len_total_rev = len (GW.reverseallPaths) 
    for path in GW.reverseallPaths:
        i_count += 1
        i_count_2 += 1
        
        cycle_inf_obj = cycle_info_optimized(path, node_dict, edge_dict, cycle_info_list)
        if cycle_inf_obj == None:
            print('An error occured')
            continue
        
        if cycle_inf_obj == 'Pass':
            #print('Exact Cycle is existed')
            continue
        
        if cycle_inf_obj.length >= min_size_of_cycle and cycle_inf_obj.length <= max_size_of_cycle and cycle_inf_obj.component_number >= min_component_number and cycle_inf_obj.component_number <= max_component_number:
            name_it += 1
            cycle_inf_obj.name = name_prefix_cycle + '_reversedOriented_' + str(name_it)
            cycle_inf_obj.reverseOriented = True
            cycle_info_list.append(cycle_inf_obj)

        print_progress_bar(i_count_2, len_total_rev, prefix='Processing:', suffix='Complete')
    
    print(cycle_info_list)

    final_str_info = ''
    print('\nCycle is reducing')
    cycle_clear_list = []

    name_it = 0
    the_other_iterator = 1
    print(len(cycle_info_list))
    
    cycle_clear_list = filter_cycles_with_kmer(cycle_info_list, k_mer_sim, threshold_sim, name_prefix_cycle)

    """
    for cycle_el in cycle_info_list:
        reverse_ori = ''
        if cycle_el.reverseOriented:
            reverse_ori = 'reverseoriented_'
        sim_num = False
        for second_seq in range(the_other_iterator,len(cycle_info_list)):
            print('.',end='',flush=True)
            if is_similar(cycle_el.sequence, cycle_info_list[second_seq].sequence, k_mer_sim, threshold_sim):
                sim_num = True
                
                break
        if sim_num == False:
            name_it += 1
            cycle_el.name = name_prefix_cycle + '_' + reverse_ori + str(name_it)
            cycle_clear_list.append(cycle_el)
        print('.',end='',flush=True)
        the_other_iterator += 1
    """
    print('\n')
    
    for cycle_el in cycle_clear_list:
        name_it += 1
        final_str_info += ">" + cycle_el.name + "-len" + str(cycle_el.length) + "-\n" + cycle_el.sequence + "\n" 
        print(cycle_el.name, 'seq_len:', cycle_el.length, ' components:', cycle_el.component_number, cycle_el.path)
    
    f = open(out_cycle_file, "w")
    f.write(final_str_info)
    f.close()

















