"""Additional taxonomy fields for genomes

Revision ID: 0a1c81836e60
Revises: 40c711d276f0
Create Date: 2016-05-26 00:28:55.264108

"""

# revision identifiers, used by Alembic.
revision = '0a1c81836e60'
down_revision = '40c711d276f0'
branch_labels = None
depends_on = None

from alembic import op
import sqlalchemy as sa


def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('genomes', schema=None) as batch_op:
        batch_op.add_column(sa.Column('tax_genus', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('tax_species', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('tax_strain', sa.String(), nullable=True))

    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('genomes', schema=None) as batch_op:
        batch_op.drop_column('tax_strain')
        batch_op.drop_column('tax_species')
        batch_op.drop_column('tax_genus')

    ### end Alembic commands ###
