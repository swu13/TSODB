import synapseclient
import pandas as pd
import os

# =========================
# 1. token
# =========================
syn = synapseclient.Synapse()
syn.login(authToken="eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTc3NzcxNzE5NiwiaWF0IjoxNzc3NzE3MTk2LCJqdGkiOiIzNjcyMiIsInN1YiI6IjMzODE1NDYifQ.d_LHkwf8nBPB-Dbmup8axmNrfaqEU27PKRH9mNFIdpgvM29rXBDGo-pctxCQA9NZJKqIF72KPkebOyjR_Njv6VlceZBfsoMYsPCFHWNGnkYwjPuL4Kdbfm-ydH8nTZdJQQhndJjdcq5O7kdvfo2HNwUHU_tY_LVKTtdNERgKlZT67Egd0xEL2dqlJirrwTab2vLr91nHShkZHyx51YCCBRHLnWXLPOsKaOeDZtyY648hJqVk-Xyvey2wIyxSWGQq6Azg0HKaKC0I1tQ5KcQ3_ec9WhX1EyZGrig4y72XwbDIZiVeN24NJqLVk3ZiH0nPuX12apeBrFVT6CHPbl465w")   

# =========================
# 2. read synapse ids
# =========================
df = pd.read_csv("synapseids.txt")   # 
ids = df["synapse_id"].dropna().unique()

# =========================
# 3. download path
# =========================
download_dir = "downloads"
os.makedirs(download_dir, exist_ok=True)

# =========================
# 4. downloading
# =========================
for sid in ids:
    try:
        print(f"Downloading {sid} ...")
        entity = syn.get(sid, downloadLocation=download_dir)
        print(f"Saved: {entity.path}")
    except Exception as e:
        print(f"? Failed {sid}: {e}")

print("? Done")