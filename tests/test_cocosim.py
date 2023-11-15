"""pytest moduel to test cocosim
"""
import pytest
import logging, os, re, subprocess
from unittest.mock import patch



#setup logger                                      
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






if __name__ == '__main__':

    pytest.main()