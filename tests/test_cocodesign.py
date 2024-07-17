import pytest
import logging
import sys
from unittest.mock import patch
from io import StringIO
from cocopaths.cocodesign import (
    set_verbosity,
    verify_domain_foldingpath,
    afp_to_domainfp,
    domain_path_to_nt_path,
    rna_design,
    extend_domain_seq,
    main
)

# Fixture to setup logger
@pytest.fixture(scope="function")
def configure_logger():
    logger = logging.getLogger("copaths")  # Use the same logger name as in your main code
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter("# %(levelname)s - %(message)s")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    logger.setLevel(logging.DEBUG)  # Set the desired logging level (e.g., DEBUG, INFO, etc.)
    yield
    logger.handlers = []

# Example content for the .pil file
PIL_CONTENT_1 = """
length a = 4
length b = 4
length c = 4
length d = 4
length e = 4
length f = 4
length g = 4
length h = 4
length L0 = 8
length S0 = 0
length S1 = 3
length S2 = 6
length S3 = 9
length S4 = 12


test_1 = L0*  S0 a L0 b S1 c* b* L0* a* d* S2 d a L0 b c S3
"""

# Helper function to create a temporary .pil file
def create_temp_pil_file(tmp_path, content):
    pil_file = tmp_path / "input.pil"
    pil_file.write_text(content)
    return pil_file

def test_main(monkeypatch, tmp_path):
    # Step 1: Create a temporary .pil file
    temp_pil_file = create_temp_pil_file(tmp_path, PIL_CONTENT_1)

    # Step 2: Mock command-line arguments
    test_args = ["prog", "-i", str(temp_pil_file), "-s", "10","--aCFP",".,(),.(),()()"]
    monkeypatch.setattr(sys, 'argv', test_args)

    # Step 3: Mock input function to simulate user input
    def mock_input(prompt=None):
        if prompt and "Please input a domain level sequence:" in prompt:
            return "L1 S1"
        elif prompt and "Please input the afp by which the domain level sequence was designed:" in prompt:
            return ".\n()\n()\n.()\n()()\n$"
        else:
            return ""  # Handle the case where prompt is None or unexpected prompt string

    monkeypatch.setattr('builtins.input', mock_input)

    print("going to main")
    nt_seq = main()

    print(f"Captured output:\n{ nt_seq = }")
    # Step 4: Verify the output and side effects
    assert nt_seq[0:8] == nt_seq[35:43]
    assert len(nt_seq) == 90


# Example content for the .pil file
PIL_CONTENT_2 = """
length a = 4
length b = 4
length c = 4
length d = 4
length e = 4
length f = 4
length g = 4
length h = 4
length L0 = 8
length S0 = 0
length S1 = 3
length S2 = 6
length S3 = 9
length S4 = 12


test_1 = a* b* L0* c* d* S0 c L0 b S1  L0*  S2 d c L0 b a S3
"""

# Helper function to create a temporary .pil file
def create_temp_pil_file(tmp_path, content):
    pil_file = tmp_path / "input.pil"
    pil_file.write_text(content)
    return pil_file

def test_main_sim_fail(monkeypatch, tmp_path):
    # Step 1: Create a temporary .pil file
    temp_pil_file = create_temp_pil_file(tmp_path, PIL_CONTENT_2)

    # Step 2: Mock command-line arguments
    test_args = ["prog", "-i", str(temp_pil_file), "-s", "10","--aCFP",".,(),.(),()()"]
    monkeypatch.setattr(sys, 'argv', test_args)

    # Step 3: Mock input function to simulate user input
    def mock_input(prompt=None):
        if prompt and "Please input a domain level sequence:" in prompt:
            return "L1 S1"
        elif prompt and "Please input the afp by which the domain level sequence was designed:" in prompt:
            return ".\n()\n()\n.()\n()()\n$"
        else:
            return ""  # Handle the case where prompt is None or unexpected prompt string

    monkeypatch.setattr('builtins.input', mock_input)

    with pytest.raises(SystemExit):
        nt_seq = main()

    
    

if __name__ == '__main__':
    pytest.main()
