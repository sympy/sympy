import csv
import os

# Path to your project root (or use current directory)
project_root = "."

# Read link replacements from CSV
replacements = {}
with open("link_replacements.csv", newline="", encoding="utf-8") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        replacements[row["Original Link"]] = row["Replacement / Updated Link"]

# Walk through all files in the project
for subdir, _, files in os.walk(project_root):
    for file in files:
        # Only modify text-based files (python, md, rst, txt, html, etc.)
        if file.endswith((".py", ".md", ".rst", ".txt", ".html")):
            file_path = os.path.join(subdir, file)
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()
            # Replace links
            modified = False
            for old, new in replacements.items():
                if old in content:
                    content = content.replace(old, new)
                    modified = True
            # Save changes if any
            if modified:
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(content)
                print(f"Updated links in {file_path}")
