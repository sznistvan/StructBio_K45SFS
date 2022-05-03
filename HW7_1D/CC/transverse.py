import matplotlib.pyplot as plt

f = open("CCHMM_cc_result.txt","r")


residues = []
prediction = []

for r in f:
	if r.startswith("SEQ"):
		print("Sequence")
		print(r)
		for i in range(5,len(r)):
			print(r[i])
			residues.append(r[i])
	if r.startswith("PRD:"):
		for j in range(5,len(r)):
			if(r[j] == 'H'):
				prediction.append(1)
			else:
				prediction.append(0)
	print("New line")

print(residues)
print(prediction)

plt.plot(prediction)
plt.show()
