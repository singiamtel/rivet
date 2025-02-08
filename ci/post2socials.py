from atproto import Client, models
from mastodon import Mastodon
import os

# Get credentials from GitLab CI/CD variables
BSKY_ID = os.getenv("BSKY_ID")
BSKY_TOKEN = os.getenv("BSKY_TOKEN")
MAST_BASE = os.getenv("MASTODON_BASE")
MAST_TOKEN = os.getenv("MASTODON_TOKEN")

# Construct the post content
MR_TITLE = os.getenv("MR_TITLE")
MR_URL = os.getenv("MR_URL")
CONTENT = f"New [Rivet analysis] - {MR_TITLE} 🎉 #HEP #OpenScience"

# Find the position of 'Rivet analysis' in the post content
start_idx = post_content.find('[Rivet analysis]')
end_idx = start_idx + len('Rivet analysis')

# Create a facet (rich text link) for the hyperlink
facet = models.AppBskyRichtextFacet.Main(
  features=[
    models.AppBskyRichtextFacet.Link(uri=MR_URL)
  ],
  index=models.AppBskyRichtextFacet.ByteSlice(byteStart=start_idx, byteEnd=end_idx)
)

# Post to Bluesky
client = Client()
client.login(BSKY_ID, BSKY_TOKEN)
client.post(CONTENT.replace('[R', 'R').replace('s]', 's'), facets=[facet])

# Post to Mastodon (doesn't support facets yet)
mastodon = Mastodon(access_token=MAST_TOKEN, api_base_url=MAST_BASE)
mastodon.status_post(CONTENT.replace('[R', 'R').replace('s]', 's'))

# Log the content
print("Submitted to socials:", CONTENT)
