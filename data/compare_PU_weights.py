import matplotlib
matplotlib.use('Agg')  # Use the Agg backend to avoid Tkinter display error
import json
import matplotlib.pyplot as plt

# Load the JSON file
with open('POG/LUM/2016preVFP_UL/puWeights.json', 'r') as f_2015:
    data_2015 = json.load(f_2015)
with open('POG/LUM/2016postVFP_UL/puWeights.json', 'r') as f_2016:
    data_2016 = json.load(f_2016)
with open('POG/LUM/2017_UL/puWeights.json', 'r') as f_2017:
    data_2017 = json.load(f_2017)
with open('POG/LUM/2018_UL/puWeights.json', 'r') as f_2018:
    data_2018 = json.load(f_2018)


nominal_weights_2015 = data_2015['corrections'][0]['data']['content'][0]['value']['content']
bin_edges_2015 = data_2015['corrections'][0]['data']['content'][0]['value']['edges']

nominal_weights_2016 = data_2016['corrections'][0]['data']['content'][0]['value']['content']
bin_edges_2016 = data_2016['corrections'][0]['data']['content'][0]['value']['edges']

nominal_weights_2017 = data_2017['corrections'][0]['data']['content'][0]['value']['content']
bin_edges_2017 = data_2017['corrections'][0]['data']['content'][0]['value']['edges']

nominal_weights_2018 = data_2018['corrections'][0]['data']['content'][0]['value']['content']
bin_edges_2018 = data_2018['corrections'][0]['data']['content'][0]['value']['edges']


# Plot the nominal and up weights
plt.figure(figsize=(10, 6))
plt.step(bin_edges_2015[:-1], nominal_weights_2015, label='2016preAPV', where='post', color='b', linestyle='-', linewidth=2)

plt.step(bin_edges_2016[:-1], nominal_weights_2016, label='2016postAPV', where='post', color='r', linestyle='-', linewidth=2)

plt.step(bin_edges_2017[:-1], nominal_weights_2017, label='2017', where='post', color='g', linestyle='-', linewidth=2)

plt.step(bin_edges_2018[:-1], nominal_weights_2018, label='2018', where='post', color='y', linestyle='-', linewidth=2)


plt.yscale('log')
#plt.step(bin_edges_2015[:-1], up_weights, label='Up', where='post', color='r', linestyle='--', linewidth=2)

# Customize the plot
plt.xlabel('Number of True Interactions')
plt.ylabel('Pileup Weight')
plt.title('Pileup Weight Distribution (2016preAPV,2016postAPV,2017,2018)')
plt.legend()
plt.grid(True)
#plt.tight_layout()

# Save the plot to a file
plt.savefig("PU_weights.png")
plt.close()  # Close the plot to free up memor