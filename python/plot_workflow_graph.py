import networkx as nx
import matplotlib.pyplot as plt

# Step 1: Parse the bash file

def parse_bash_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    commands = []
    loop_stack = []
    for line in lines:
        line = line.strip()
        if line.startswith('#') or not line:
            continue

        if line.startswith('for ') or line.startswith('while '):
            loop_stack.append(line)
            commands.append(line)
        elif line.startswith('done'):
            loop_start = loop_stack.pop()
            commands.append(f'end {loop_start}')
        else:
            commands.append(line)

    return commands

# Step 2: Build the workflow tree
def build_workflow_tree(commands):
    G = nx.DiGraph()
    prev_command = None

    for command in commands:
        if '&&' in command:
            parts = command.split('&&')
        elif '||' in command:
            parts = command.split('||')
        elif '|' in command:
            parts = command.split('|')
        else:
            parts = [command]

        for part in parts:
            part = part.strip()
            if prev_command:
                G.add_edge(prev_command, part)
            prev_command = part
    
    return G

# Step 3: Plot the workflow tree
def plot_workflow_tree(G):
    pos = nx.spring_layout(G)
    plt.figure(figsize=(10, 8))
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=2000, font_size=10, font_weight='bold', arrows=True)
    plt.title('Workflow Tree')
    plt.show()

# Main function
def main(filepath):
    commands = parse_bash_file(filepath)
    G = build_workflow_tree(commands)
    plot_workflow_tree(G)

# Example usage
filepath = 'fit_normalisation.sh'
main(filepath)
