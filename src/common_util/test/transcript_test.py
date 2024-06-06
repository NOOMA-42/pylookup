import unittest
from src.common_util.transcript import CommonTranscript
from src.common_util.curve import Scalar

class TestTranscript(unittest.TestCase):
    def test_transcript(self):
        scalar = Scalar(123)
        transcript = CommonTranscript(b"plonk")
        transcript.append_scalar(b"scalar", scalar)
        transcript2 = CommonTranscript(b"plonk")
        transcript2.append_scalar(b"scalar", scalar)
        self.assertEqual(
            transcript.get_and_append_challenge(b"scalar"), 
            transcript2.get_and_append_challenge(b"scalar")     
        )

        