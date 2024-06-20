import json
import os

import openai
from dotenv import load_dotenv

load_dotenv()

openai.api_key = os.getenv("OPENAI_API_KEY")

def clean_article_text(input_text):
    # Count the number of tokens in the input text
    prompt_tokens = len(input_text.split())  # This is a simple approximation

    # Set max_tokens for the output to be roughly the same as the input or slightly less
    max_response_tokens = max(1, prompt_tokens)

    # Define the prompt to clean the article text
    prompt = f"""
    Extract the main content of the following text and remove any unnecessary information such as login information, subscription details, cookies notices, etc. Keep only the essential parts of the article:

    {input_text}

    Cleaned Article Text:
    """

    # Call the OpenAI API with the prompt
    response = openai.Completion.create(
        model="text-davinci-003",  # or another model of your choice
        prompt=prompt,
        temperature=0,
        max_tokens=max_response_tokens,  # Adjust as necessary
        top_p=1,
        frequency_penalty=0,
        presence_penalty=0
    )

    # Extract and return the cleaned text from the response
    cleaned_text = response.choices[0].text.strip()
    return cleaned_text

def main():
    try:
        with open('extracted_text_1.txt', 'r', encoding='utf-8') as file:
            input_text = file.read()

        # Get the cleaned article text
        cleaned_text = clean_article_text(input_text)

        # Print the cleaned text
        print(cleaned_text)
    except Exception as e:
        print(f"An error occurred: {e}")
        
if __name__ == "__main__":
    main()