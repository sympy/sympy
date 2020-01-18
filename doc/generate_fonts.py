import os

sympy_root_dir = os.path.abspath(os.path.join(__file__, '..', '..'))
sympy_dir = os.path.join(sympy_root_dir, 'sympy')
sympy_doc_dir = os.path.join(sympy_root_dir, 'doc')

font_path = os.path.join(sympy_doc_dir, 'fonts', "DejaVuSansMono.ttf")
output_dir = os.path.join(sympy_doc_dir, '_build', 'html', '_static')
output_font_path = os.path.join(output_dir, "DejaVuSansMono.woff")


def get_char_set_from_file(filename: str):
    char_set = set()

    with open(filename, 'r', encoding="utf-8") as f:
        new = set(f.read())
        char_set = char_set.union(new)

    return char_set

def get_char_set_from_sympy():
    char_set = set()

    for path, folders, files in os.walk(sympy_dir):
        for file in files:
            filename, ext = os.path.splitext(file)
            if ext == '.py':
                new = get_char_set_from_file(os.path.join(path, filename+ext))
                char_set = char_set.union(new)

    for path, folders, files in os.walk(sympy_doc_dir):
        for file in files:
            filename, ext = os.path.splitext(file)
            if ext == '.rst':
                new = get_char_set_from_file(os.path.join(path, filename+ext))
                char_set = char_set.union(new)

    return char_set


if __name__ == '__main__':
    char_set = get_char_set_from_sympy()
    text = "".join(char_set)

    import subprocess
    subprocess.run(
        [
            'pyftsubset', font_path,
            '--text="{}"'.format(text),
            '--output-file={}'.format(output_font_path)
        ]
    )
