# Stage 1: crop a local H&E context patch per Visium spot and embed it with the
# Phikon pathology foundation model (Owkin, iBOT/ViT-B trained on TCGA histology).
# NB: only the hires thumbnail (~8.5 px/spot) is available, so we use a context
# window (~9-spot neighbourhood) — a resolution-limited proof of concept.
import os, json, numpy as np
from PIL import Image
import torch
from transformers import AutoImageProcessor, AutoModel
ROOT="/scratch/chahal.d/Spatial_GenAI"; OUT=os.path.join(ROOT,"HE_FoundationModel")
SAMPLES=["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]
WIN=48   # half-window in hires pixels (~9-spot neighbourhood)
dev="cuda" if torch.cuda.is_available() else "cpu"; print("device",dev, flush=True)
proc=AutoImageProcessor.from_pretrained("owkin/phikon")
model=AutoModel.from_pretrained("owkin/phikon").to(dev).eval()

embs=[]; meta=[]
for s in SAMPLES:
    sp=os.path.join(ROOT,"samples",s,"spatial")
    scal=json.load(open(os.path.join(sp,"scalefactors_json.json")))["tissue_hires_scalef"]
    img=Image.open(os.path.join(sp,"tissue_hires_image.png")).convert("RGB"); W,H=img.size
    rows=[l.strip().split(",") for l in open(os.path.join(sp,"tissue_positions_list.csv"))]
    patches=[]; bcs=[]
    for r in rows:
        if r[1]!="1": continue
        yc=float(r[4])*scal; xc=float(r[5])*scal
        x0=max(0,int(xc-WIN)); y0=max(0,int(yc-WIN)); x1=min(W,int(xc+WIN)); y1=min(H,int(yc+WIN))
        patches.append(img.crop((x0,y0,x1,y1))); bcs.append(r[0])
    for i in range(0,len(patches),128):
        inp=proc(images=patches[i:i+128], return_tensors="pt").to(dev)
        with torch.no_grad(): out=model(**inp)
        embs.append(out.last_hidden_state[:,0,:].cpu().numpy())
    meta += [(s,b) for b in bcs]; print(s, len(bcs), "spots", flush=True)
X=np.concatenate(embs,0); meta=np.array(meta)
np.savez(os.path.join(OUT,"phikon_embeddings.npz"), X=X, sample=meta[:,0], barcode=meta[:,1])
print("EMB_DONE", X.shape, flush=True)
