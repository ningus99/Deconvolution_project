
import numpy as np
import pandas as pd 
ruta_archivo = "/storage/kuijjerarea/mario/sisana/network/Predicted/norm/results/lioness.npy"


data = np.load(ruta_archivo)

df = pd.DataFrame(data)

print("Type of data:", type(data))


if isinstance(data, np.ndarray):
    print("Array_shape:", data.shape)  
    print("Element_type:", data.dtype)  
    print("First_elements:\n", data[:5])  

df.to_csv("raw_edges_pred.csv", index=True)
