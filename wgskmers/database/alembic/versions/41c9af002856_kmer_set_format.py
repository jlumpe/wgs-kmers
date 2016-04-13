"""Kmer set format

Revision ID: 41c9af002856
Revises: 704356629cab
Create Date: 2016-04-18 21:49:06.447250

"""

# revision identifiers, used by Alembic.
revision = '41c9af002856'
down_revision = '704356629cab'
branch_labels = None
depends_on = None

from alembic import op
import sqlalchemy as sa


def upgrade():
    op.add_column('kmer_collections', sa.Column('format', sa.String(), nullable=True))
    op.execute("UPDATE kmer_collections SET format = 'raw'")
    with op.batch_alter_table('kmer_collections', schema=None) as batch_op:
        batch_op.alter_column('format', nullable=False)


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('kmer_collections', schema=None) as batch_op:
        batch_op.drop_column('format')

    ### end Alembic commands ###