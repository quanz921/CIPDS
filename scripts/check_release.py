"""Verify the public release boundary without opening external research inputs."""
from pathlib import Path
import ast,csv,hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[1]
allowed_suffixes={".R",".py",".md",".txt",".json",".csv"}
allowed_csv={"demo/synthetic_preview.csv","demo/expected_demo_metrics.csv"}
fail=[];files=[]
secret_patterns=[
 re.compile(r"\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{20,}\b"),
 re.compile(r"\bgithub_pat_[A-Za-z0-9_]{30,}\b"),
 re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
 re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
 re.compile(r"(?i)\b(?:api_key|password|access_token)\s*[:=]\s*['\"][^'\"\s]{12,}['\"]"),
]
for f in ROOT.rglob("*"):
 if not f.is_file():continue
 rel=f.relative_to(ROOT).as_posix()
 if rel.startswith("demo_work/") or "__pycache__" in f.parts or ".git" in f.relative_to(ROOT).parts:continue
 files.append(rel)
 if f.is_symlink():fail.append(rel+": symlink forbidden")
 if f.name not in {".gitignore",".gitattributes"} and f.suffix not in allowed_suffixes:fail.append(rel+": unapproved file type")
 if f.stat().st_size>2_000_000:fail.append(rel+": unexpectedly large")
 try:text=f.read_text(encoding="utf-8-sig")
 except UnicodeError:fail.append(rel+": binary/non-UTF8 content");continue
 if f.suffix==".py":
  try:ast.parse(text)
  except SyntaxError as e:fail.append(rel+": Python syntax error "+str(e.lineno))
 if f.suffix==".csv":
  if rel not in allowed_csv:fail.append(rel+": CSV not explicitly allowed")
  else:
   with f.open(encoding="utf-8-sig",newline="") as h:rows=list(csv.DictReader(h))
   if not rows or any(r.get("is_synthetic","").lower() not in ("true","1") for r in rows):
    fail.append(rel+": missing per-row synthetic marker")
   if "preview" in rel and any(not r.get("synthetic_id","").startswith("DEMO_") for r in rows):fail.append(rel+": unexpected identifier")
 if rel!="scripts/check_release.py":
  if any(p.search(text) for p in secret_patterns):fail.append(rel+": possible credential literal (content withheld)")
  if re.search(r"[CD]:[/\\](?:Users|CodexData)",text,re.I):fail.append(rel+": author-machine absolute path")
manifest=json.loads((ROOT/"provenance/source_manifest.json").read_text(encoding="utf-8-sig"))
for item in manifest:
 path=ROOT/item["file"]
 if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=item["published_sha256"]:fail.append(item["file"]+": source checksum mismatch")
release=ROOT/"provenance/release_manifest.json"
if release.exists():
 entries=json.loads(release.read_text(encoding="utf-8-sig"))["files"]
 if set(x["file"] for x in entries)!=set(files)-{"provenance/release_manifest.json"}:fail.append("Release file set differs from manifest")
 for x in entries:
  if not (ROOT/x["file"]).exists() or hashlib.sha256((ROOT/x["file"]).read_bytes()).hexdigest()!=x["sha256"]:fail.append(x["file"]+": release checksum mismatch")
print(json.dumps({"scope":"public files; demo_work is excluded","files":len(files),"source_files":len(manifest),"failures":fail,"passed":not fail},indent=2))
sys.exit(bool(fail))

