import os
import markdown

def convert_markdown_to_html(markdown_file_path, output_html_path):
    # Read the Markdown file
    with open(markdown_file_path, 'r', encoding='utf-8') as f:
        md_content = f.read()

    # Convert Markdown to HTML
    html_content = markdown.markdown(md_content)

    # Save the HTML content to an output file
    with open(output_html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

# Define paths

script_dir = os.path.dirname(os.path.abspath(__file__))
folder_case="autotest/"
folder_case = os.path.join(script_dir, "../"+folder_case)

markdown_file = folder_case+"autotest.md"
output_html = folder_case+"autotest.html"

# Convert the Markdown file to HTML
convert_markdown_to_html(markdown_file, output_html)